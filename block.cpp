/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Sebastian Deorowicz and Szymon Grabowski
  
  Version: 0.2
*/

#include "block.h"
#include "superblock.h"
#include <algorithm>

using namespace std;

const uint32 CBlock::MAX_SIZE = 32;

// ********************************************************************************************
//
// ********************************************************************************************
CBlock::CBlock(uint64 _b_start_rec_no, uint32 _rec_count) 
{
	b_start_rec_no = _b_start_rec_no;
	rec_count      = _rec_count;

	qua_stream = NULL;
	run_stream = NULL;
	
	records.resize(MAX_SIZE);
};

// ********************************************************************************************
CBlock::~CBlock()
{
	delete[] qua_stream;
	delete[] run_stream;
}

// ********************************************************************************************
void CBlock::Reset(uint64 _b_start_rec_no, uint32 _global_max_seq_length, uint32 _rec_count)
{
	b_start_rec_no        = _b_start_rec_no;
	rec_count             = _rec_count;
	global_max_seq_length = _global_max_seq_length;

	delete[] qua_stream;
	delete[] run_stream;
	qua_stream = NULL;
	run_stream = NULL;
}

// ********************************************************************************************
void CBlock::InsertRecord(const CFastqRecord &rec, vector<uint32> &dna_occ, bool &is_plus_only, vector<uint32> &quality_occ)
{
	for(uint32 i = 0; i < rec.seq_len; ++i)
		++dna_occ[rec.seq[i]];

	is_plus_only &= rec.plus_len == 1;

	records[rec_count++].Set(rec);
}

// ********************************************************************************************
bool CBlock::Process(CBitStream &bit_stream, CLZMatcher &lz_matcher, vector<CField> &fields, uint32 n_fields,
	bool _is_plus_only, bool _is_dna_plain, bool _is_length_variable,
	unsigned char *sym_code, unsigned char *qua_code, CHuffman::t_code **qua_huf_codes, 
	uint32 max_run_len, CHuffman::t_code **run_huf_codes, uint32 n_qualities, uint32 block_no, 
	uint32 _quality_stats_mode, uint32 _try_lz)
{
	is_dna_plain	   = _is_dna_plain;
	quality_stats_mode = _quality_stats_mode;
	try_lz			   = _try_lz;

	uint32 i;
	if(!_is_plus_only)
		for(i = 0; i < rec_count; ++i)
			bit_stream.PutBit(records[i].plus_len == 1);
	uint32 quality_len_bits = bit_stream.BitLength(global_max_seq_length);
	if(_is_length_variable)
		for(i = 0; i < rec_count; ++i)
			bit_stream.PutBits(records[i].quality_len, quality_len_bits);
	bit_stream.FlushPartialWordBuffer();

	StoreTitle(bit_stream, fields, block_no);

	if(quality_stats_mode == 2)
		StoreQualityRLE(bit_stream, qua_code, qua_huf_codes, run_huf_codes, n_qualities);
	else
		StoreQualityPlain(bit_stream, qua_code, qua_huf_codes, n_qualities);

	StoreDNAPlain(bit_stream, lz_matcher, _is_dna_plain, sym_code);

	return true;
}

// ********************************************************************************************
bool CBlock::Read(CBitStream &bit_stream, CLZMatcher &lz_matcher, vector<CField> &fields, uint32 n_fields,
	bool _is_plus_only, bool _is_dna_plain, bool _is_length_variable,
	vector<unsigned char> &symbols, vector<unsigned char> &qualities, vector<CHuffman*> &Huffman_qua, 
	uint32 max_run_len, vector<CHuffman*> &Huffman_run, uint32 n_qualities, uint32 block_no, 
	uint32 _quality_stats_mode, uint32 _try_lz, bool extracting)
{
	uint32 i;

	is_dna_plain       = _is_dna_plain;
	quality_stats_mode = _quality_stats_mode;
	try_lz			   = _try_lz;

	no_of_N.resize(rec_count);
	for(i = 0; i < rec_count; ++i)
	{
		records[i].Reset();
		no_of_N[i] = 0;
	}

	int32 quality_len_bits = bit_stream.BitLength(global_max_seq_length);
	if(!_is_plus_only)
		for(i = 0; i < rec_count; ++i)
		{
			bit_stream.GetBit(records[i].plus_len);
		}
	else
		for(i = 0; i < rec_count; ++i)
			records[i].plus_len = 1;
	if(_is_length_variable)
		for(i = 0; i < rec_count; ++i)
		{
			uint32 tmp;
			bit_stream.GetBits(tmp, quality_len_bits);
			records[i].quality_len = tmp;
			records[i].ExtendTo(records[i].quality, records[i].quality_size, tmp+2);
			records[i].seq_len = tmp;
			records[i].ExtendTo(records[i].seq, records[i].seq_size, tmp+2);
		}
	else
	{
		for(i = 0; i < rec_count; ++i)
		{
			records[i].quality_len = global_max_seq_length;
			records[i].ExtendTo(records[i].quality, records[i].quality_size, global_max_seq_length+2);
			records[i].seq_len = global_max_seq_length;
			records[i].ExtendTo(records[i].seq, records[i].seq_size, global_max_seq_length+2);
		}
	}
	bit_stream.FlushInputWordBuffer();

	ReadTitle(bit_stream, fields, n_fields, block_no);

	if(quality_stats_mode == 2)
	{
		ReadQualityRLE(bit_stream, qualities, Huffman_qua, n_qualities, Huffman_run, max_run_len);
		MakeUnRLE();
	}
	else
		ReadQualityPlain(bit_stream, qualities, Huffman_qua, n_qualities);

	ReadDNAPlain(bit_stream, lz_matcher, true, symbols, extracting);

	return true;
}

// ********************************************************************************************
bool CBlock::StoreTitle(CBitStream &bit_stream, vector<CField> &fields, int32 block_no)
{
	uint32 i;
	uint32 n_fields = (uint32) fields.size();
	prev_value.resize(n_fields);

	for(i = 0; i < n_fields; ++i)
	{
		if(fields[i].is_constant)
			continue;
		prev_value[i] = 0;
		if(!fields[i].is_numeric)
			bit_stream.PutBit(fields[i].block_desc[block_no].is_block_constant);
		if(fields[i].is_numeric)
		{
			fields[i].block_desc[block_no].is_block_delta_constant &= fields[i].block_desc[block_no].block_delta_constant == fields[i].min_delta;
			if(fields[i].is_delta_coding)
				bit_stream.PutBit(fields[i].block_desc[block_no].is_block_delta_constant);
			else
				bit_stream.PutBit(fields[i].block_desc[block_no].is_block_value_constant);
		}		
	}

	for(i = 0; i < rec_count; ++i)
	{
		uint32 c_field = 0;
		uint32 start_pos = 0;
		CFastqRecord &rec = records[i];
		for(uint32 k = 0; k <= rec.title_len; ++k)
		{
			if(rec.title[k] != fields[c_field].sep)
				continue;

			CField &cur_field = fields[c_field];

			if(cur_field.is_constant)
			{
				start_pos = k+1;
				c_field++;
				continue;
			}
			if(cur_field.is_numeric)
			{
				int32 value = cur_field.to_num(rec.title+start_pos, k-start_pos);
				if(i == 0)
					bit_stream.PutBits(value-cur_field.min_value, cur_field.no_of_bits_per_value);
				else if((cur_field.is_delta_coding && !cur_field.block_desc[block_no].is_block_delta_constant) ||
					(!cur_field.is_delta_coding && !cur_field.block_desc[block_no].is_block_value_constant))
				{
					int32 to_store;
					if(cur_field.is_delta_coding)
						to_store = value - prev_value[c_field] - cur_field.min_delta;
					else
						to_store = value - fields[c_field].min_value;
					if(cur_field.global_stats)
						bit_stream.PutBits(cur_field.Huffman_global->codes[to_store].code, cur_field.Huffman_global->codes[to_store].len);
					else
						bit_stream.PutBits(to_store, cur_field.no_of_bits_per_num);
				}

				prev_value[c_field] = value;
				start_pos = k+1;
				c_field++;
				continue;
			}
			if(i > 0 && cur_field.block_desc[block_no].is_block_constant)
			{
				start_pos = k+1;
				c_field++;
				continue;
			}
			if(!cur_field.is_len_constant)
				bit_stream.PutBits(k-start_pos - cur_field.min_len, cur_field.no_of_bits_per_len);
			for(uint32 j = 0; j < k-start_pos; ++j)
				if(j < cur_field.len && cur_field.Ham_mask[j])
					;
				else
				{
					unsigned char c = rec.title[start_pos+j];
					bit_stream.PutBits(cur_field.Huffman_local[MIN(j, CSuperBlock::MAX_FIELD_STAT_LEN)]->codes[c].code, 
						cur_field.Huffman_local[MIN(j, CSuperBlock::MAX_FIELD_STAT_LEN)]->codes[c].len);
				}

			start_pos = k+1;
			c_field++;
		}
	}
	bit_stream.FlushPartialWordBuffer();

	return true;
}

// ********************************************************************************************
bool CBlock::StoreQualityPlain(CBitStream &bit_stream, unsigned char *qua_code, CHuffman::t_code **qua_huf_codes, int32 n_qualities)
{
	// Quality data
	for(uint32 i = 0; i < rec_count; ++i)
	{
		unsigned char *cur_quality = records[i].quality; 
		uint32 cur_quality_len = records[i].quality_len;
		uint32 cur_quality_len_th;

		if(quality_stats_mode == 1)
			cur_quality_len_th = cur_quality_len;
		else
			cur_quality_len_th = records[i].rec_th_len;

		if(quality_stats_mode == 3)			// Truncate #
		{
			bit_stream.PutBits(cur_quality_len_th != cur_quality_len, 1);
			if(cur_quality_len_th != cur_quality_len)
				bit_stream.PutBits(cur_quality_len - cur_quality_len_th, bit_stream.BitLength(cur_quality_len));
		}

		for(uint32 j = 0; j < cur_quality_len_th; ++j)
		{
			int32 qua = qua_code[cur_quality[j]];
			bit_stream.PutBits(qua_huf_codes[j+1][qua].code, qua_huf_codes[j+1][qua].len);
		}
	}
	bit_stream.FlushPartialWordBuffer();

	return true;
}

// ********************************************************************************************
bool CBlock::StoreQualityRLE(CBitStream &bit_stream, unsigned char *qua_code, CHuffman::t_code **qua_huf_codes, CHuffman::t_code **run_huf_codes, 
	int32 n_qualities)
{
	// Quality data
	unsigned char prev = 0;
	uint32 pos = 0;

	for(uint32 i = 0; i < qua_stream_len; ++i)
	{
		unsigned char qua = qua_code[qua_stream[i]];
		unsigned char len = run_stream[i];
		bit_stream.PutBits(qua_huf_codes[prev][qua].code, qua_huf_codes[prev][qua].len);
		bit_stream.PutBits(run_huf_codes[qua][len].code, run_huf_codes[qua][len].len);
		prev = qua;
		pos += len+1;
	}

	bit_stream.FlushPartialWordBuffer();

	return true;
}

// ********************************************************************************************
bool CBlock::StoreDNAPlain(CBitStream &bit_stream, CLZMatcher &lz_matcher, bool _is_dna_plain, unsigned char *sym_code)
{
	uint32 i;

	// Info about LZ matches
	FindLZMatches(lz_matcher);
	if(try_lz)
	{
		for(i = 0; i < rec_count; ++i)
		{
			bit_stream.PutBit(lz_matches[i].length != 0);
			if(lz_matches[i].length == 0)
				bit_stream.PutBit(records[i].lz_inserted);
		}
		bit_stream.FlushPartialWordBuffer();
	
		// DNA data
		uint32 rec_no_bits = bit_stream.BitLength(b_start_rec_no + rec_count - 1);
		uint32 offset_bits = bit_stream.BitLength(global_max_seq_length - CLZMatcher::MIN_MATCH_LEN);
		uint32 length_bits;
		for(i = 0; i < rec_count; ++i)
		{
			if(lz_matches[i].length > 0)
			{
				length_bits = bit_stream.BitLength(MIN(records[i].seq_len - CLZMatcher::MIN_MATCH_LEN, 255));
	
				bit_stream.PutBits(lz_matches[i].rec_no, rec_no_bits);
				bit_stream.PutBits(lz_matches[i].length - CLZMatcher::MIN_MATCH_LEN, length_bits);
				bit_stream.PutBits(lz_matches[i].rec_offset, offset_bits);
			}
		}
		bit_stream.FlushPartialWordBuffer();
	}

	for(i = 0; i < rec_count; ++i)
	{
		uint32 cur_seq_len = records[i].seq_len;
		unsigned char *cur_seq = records[i].seq;
		if(is_dna_plain)
			for(uint32 j = lz_matches[i].length; j < cur_seq_len; ++j)
				bit_stream.Put2Bits(sym_code[cur_seq[j]]);
	}
	bit_stream.FlushPartialWordBuffer();

	return true;
}

// ********************************************************************************************
bool CBlock::ReadTitle(CBitStream &bit_stream, vector<CField> &fields, uint32 n_fields, int32 block_no)
{
	uint32 i;
	prev_value.resize(n_fields);

	for(i = 0; i < n_fields; ++i)
	{
		uint32 tmp;
		if(fields[i].is_constant)
			continue;
		prev_value[i] = 0;
		if(!fields[i].is_numeric)
		{
			bit_stream.GetBit(tmp);
			fields[i].block_desc[block_no].is_block_constant = tmp != 0;
		}
		if(fields[i].is_numeric)
		{
			bit_stream.GetBit(tmp);
			if(fields[i].is_delta_coding)
				fields[i].block_desc[block_no].is_block_delta_constant = tmp != 0;
			else
				fields[i].block_desc[block_no].is_block_value_constant = tmp != 0;
		}		
	}

	for(i = 0; i < rec_count; ++i)
	{
		uint32 c_field = 0;
		uint32 start_pos = 0;

		for(uint32 j = 0; j < n_fields; ++j)
		{
			CField &cur_field = fields[j];
			if(cur_field.is_constant)
			{
				records[i].AppendTitle(cur_field.data, cur_field.len);
				records[i].AppendTitle(cur_field.sep);
				continue;
			}
			if(cur_field.is_numeric)
			{
				uint32 tmp;
				if(records[i].title_len + 10 >= records[i].title_size)
					records[i].Extend(records[i].title, records[i].title_size);
				if(i == 0)
				{
					bit_stream.GetBits(tmp, cur_field.no_of_bits_per_value);
					tmp += cur_field.min_value;
					records[i].title_len += cur_field.to_string(records[i].title+records[i].title_len, tmp);
					prev_value[j] = tmp;
				}
				else
				{
					if((cur_field.is_delta_coding && !cur_field.block_desc[block_no].is_block_delta_constant) ||
						(!cur_field.is_delta_coding && !cur_field.block_desc[block_no].is_block_value_constant))
					{
						if(cur_field.no_of_bits_per_num)
						{
							if(cur_field.global_stats)
							{
								uint32 bit;
								int32 h_tmp;
								do
								{
									bit_stream.GetBit(bit);
									h_tmp = cur_field.Huffman_global->Decode(bit);
								} while(h_tmp < 0);
								tmp = h_tmp;
							}
							else
								bit_stream.GetBits(tmp, cur_field.no_of_bits_per_num);
						}
						else
							tmp = 0;
					}
					else
					{
						if(cur_field.is_delta_coding)
							tmp = 0;
						else
							tmp = prev_value[j] - cur_field.min_value;
					}
					if(cur_field.is_delta_coding)
						tmp += prev_value[j] + cur_field.min_delta;
					else
						tmp += cur_field.min_value;
					records[i].title_len += cur_field.to_string(records[i].title+records[i].title_len, tmp);
					prev_value[j] = tmp;
				}
				records[i].AppendTitle(cur_field.sep);
				continue;
			}

			if(i > 0 && cur_field.block_desc[block_no].is_block_constant)
			{
				records[i].AppendTitle(records[0].title+cur_field.block_str_start, cur_field.block_str_len);
				records[i].AppendTitle(cur_field.sep);
				continue;
			}

			uint32 field_len;
			if(!cur_field.is_len_constant)
			{
				bit_stream.GetBits(field_len, cur_field.no_of_bits_per_len);			
				field_len += cur_field.min_len;
			}
			else
				field_len = cur_field.len;
			
			for(uint32 k = 0; k < field_len; ++k)
			{
				if(k < cur_field.len && cur_field.Ham_mask[k])
					records[i].AppendTitle(cur_field.data[k]);
				else
				{
					uint32 bit;
					uint32 tmp;
					int32 h_tmp;
					CHuffman *cur_huf = cur_field.Huffman_local[MIN(k, CSuperBlock::MAX_FIELD_STAT_LEN)]; 
					do
					{
						bit_stream.GetBit(bit);
						h_tmp = cur_huf->Decode(bit);
					} while(h_tmp < 0);
					tmp = h_tmp;
					records[i].AppendTitle((unsigned char) tmp);
				}
			}
			if(i == 0 && cur_field.block_desc[block_no].is_block_constant)
			{
				cur_field.block_str_start = records[i].title_len - field_len;
				cur_field.block_str_len   = field_len;
			}

			records[i].AppendTitle(cur_field.sep);
		}

		records[i].title_len--;			// do not count last '\0' symbols
		records[i].plus[0] = '+';
		if(records[i].plus_len)
			records[i].plus[1] = '\0';
		else
		{
			records[i].ExtendTo(records[i].plus, records[i].plus_len, records[i].title_len+1);
			copy(records[i].title+1, records[i].title+records[i].title_len, records[i].plus+1);
		}
	}

	bit_stream.FlushInputWordBuffer();
	return true;
}

// ********************************************************************************************
bool CBlock::ReadQualityPlain(CBitStream &bit_stream, vector<unsigned char> &qualities, 
	vector<CHuffman*> &Huffman_qua, int32 n_qualities)
{
	// Quality data
	for(uint32 i = 0; i < rec_count; ++i)
	{
		no_of_N[i] = 0;
		unsigned char *cur_quality = records[i].quality; 
		uint32 cur_quality_len = records[i].quality_len;
		uint32 trunc_len = 0;

		if(quality_stats_mode == 3)			// Truncate #
		{
			uint32 is_truncated;
			bit_stream.GetBit(is_truncated);
			if(is_truncated)
				bit_stream.GetBits(trunc_len, bit_stream.BitLength(cur_quality_len));
		}

		for(uint32 j = 0; j < cur_quality_len-trunc_len; ++j)
		{
			uint32 bit;
			int32 h_tmp;
			do
			{
				bit_stream.GetBit(bit);
				h_tmp = Huffman_qua[j+1]->Decode(bit);
			} while(h_tmp < 0);

			if((cur_quality[j] = qualities[h_tmp]) >= 128)
				no_of_N[i]++;
		}
		for(uint32 j = cur_quality_len-trunc_len; j < cur_quality_len; ++j)
			cur_quality[j] = '#';
	}
	bit_stream.FlushInputWordBuffer();

	return true;
}

// ********************************************************************************************
bool CBlock::ReadQualityRLE(CBitStream &bit_stream, vector<unsigned char> &qualities, 
vector<CHuffman*> &Huffman_qua, int32 n_qualities, vector<CHuffman*> &Huffman_run, int32 max_run_len)
{
	uint32 i;
	int32 max_len = 0;
	for(i = 0; i < rec_count; ++i)
		max_len += records[i].quality_len;

	qua_stream_len = run_stream_len = max_len;
	qua_stream = new unsigned char[max_len];
	run_stream = new unsigned char[max_len];
	int32 qua_stream_idx = 0;
	int32 run_stream_idx = 0;

	uint32 prev = 0;
	i = 0;
	for(uint32 j = 0; j < qua_stream_len; ++i)
	{
		int32 h_tmp;
		uint32 bit;
		// Quality data
		do
		{
			bit_stream.GetBit(bit);
			h_tmp = Huffman_qua[prev]->Decode(bit);
		} while (h_tmp < 0);
		qua_stream[i] = qualities[h_tmp];
		prev = h_tmp;

		// Run data
		do
		{
			bit_stream.GetBit(bit);
			h_tmp = Huffman_run[prev]->Decode(bit);
		} while(h_tmp < 0);
		run_stream[i] = (unsigned char) h_tmp;

		j += h_tmp+1;
	}

	bit_stream.FlushInputWordBuffer();

	return true;
}

// ********************************************************************************************
bool CBlock::ReadDNAPlain(CBitStream &bit_stream, CLZMatcher &lz_matcher, bool _is_dna_plain, 
	vector<unsigned char> &symbols, bool extracting)
{
	uint32 i;

	// Info about LZ matches
	lz_matches.resize(rec_count);
	
	for(i = 0; i < rec_count; ++i)
		records[i].seq_len = records[i].quality_len - no_of_N[i];

	if(try_lz)
	{
		for(i = 0; i < rec_count; ++i)
		{
			uint32 tmp;
			bit_stream.GetBit(tmp);
			lz_matches[i].length = tmp;
			if(!lz_matches[i].length)
			{
				bit_stream.GetBit(tmp);
				records[i].lz_inserted = tmp != 0;
			}
		}
		bit_stream.FlushInputWordBuffer();

		// DNA data
		uint32 rec_no_bits = bit_stream.BitLength(b_start_rec_no + rec_count - 1);
		uint32 offset_bits = bit_stream.BitLength(global_max_seq_length - CLZMatcher::MIN_MATCH_LEN);
		uint32 length_bits;

		for(i = 0; i < rec_count; ++i)
		{
			if(lz_matches[i].length > 0)
			{
				uint32 tmp;
				bit_stream.GetBits(tmp, rec_no_bits);
				lz_matches[i].rec_no = tmp;

				length_bits = bit_stream.BitLength(MIN(records[i].seq_len - CLZMatcher::MIN_MATCH_LEN, 255));
				bit_stream.GetBits(tmp, length_bits);
				lz_matches[i].length = tmp + CLZMatcher::MIN_MATCH_LEN;
				bit_stream.GetBits(tmp, offset_bits);
				lz_matches[i].rec_offset = tmp;
			}
		}
		bit_stream.FlushInputWordBuffer();
	}

	for(i = 0; i < rec_count; ++i)
	{
		if(try_lz)
		{
			if(!extracting)
				DecodeLZMatches(lz_matcher, i);
		}
		else
			lz_matches[i].length = 0;

		uint32 cur_seq_len = records[i].seq_len;
		unsigned char *cur_seq = records[i].seq;
		if(is_dna_plain)
		{
			uint32 j;
			for(j = lz_matches[i].length; j < cur_seq_len; ++j)
			{
				uint32 tmp;
				bit_stream.Get2Bits(tmp);
				cur_seq[j] = symbols[tmp];
			}	
			cur_seq[j] = '\0';
		}
		if(try_lz && !extracting)
			DecodeLZMatches_Insert(lz_matcher, i);
	}
	bit_stream.FlushInputWordBuffer();

	return true;
}

// ********************************************************************************************
void CBlock::FindLZMatches(CLZMatcher &lz_matcher)
{
	lz_matches.resize(records.size());
	for(uint32 i = 0; i < records.size(); ++i)
		lz_matches[i].length = 0;

	if(!try_lz)
		return;

	for(uint32 i = 0; i < rec_count; ++i)
	{
		uint32 rec_no;
		uint32 rec_offset;
		uint32 match_len;
		bool match_found = lz_matcher.FindMatch(records[i].seq, records[i].seq_len, records[i].quality, records[i].quality_len, 
			rec_no, rec_offset, match_len);

		if(match_found)
		{
			lz_matches[i].rec_no      = rec_no;
			lz_matches[i].rec_offset  = rec_offset;
			lz_matches[i].length      = match_len;
			records[i].lz_inserted    = false;
		}
		else
			records[i].lz_inserted = lz_matcher.InsertEncoding(b_start_rec_no+i, records[i].seq, records[i].seq_len, records[i].quality, records[i].quality_len);
	}
}

// ********************************************************************************************
void CBlock::DecodeLZMatches(CLZMatcher &lz_matcher, uint32 rec_no)
{
	if(!try_lz)
		return;
		
	if(lz_matches[rec_no].length)
		lz_matcher.DecodeMatch(records[rec_no].seq, lz_matches[rec_no].rec_no, lz_matches[rec_no].rec_offset, lz_matches[rec_no].length,
			records[rec_no].seq_len);
}

// ********************************************************************************************
void CBlock::DecodeLZMatches_Insert(CLZMatcher &lz_matcher, uint32 rec_no)
{
	if(!try_lz)
		return;

	if(records[rec_no].lz_inserted)
		lz_matcher.InsertDecoding(b_start_rec_no+rec_no, records[rec_no].seq, records[rec_no].seq_len, records[rec_no].quality, 
			records[rec_no].quality_len);
}

// ********************************************************************************************
bool CBlock::MakeRLE()
{
	uint32 i, j;
	uint32 tot_qua_len = 0;
	
	for(i = 0; i < rec_count; ++i)
		tot_qua_len += records[i].quality_len;

	qua_stream = new unsigned char[tot_qua_len];
	run_stream = new unsigned char[tot_qua_len];

	unsigned char prev = 0;
	uint32 run_len = 0;
	uint32 qua_idx = 0;
	uint32 run_idx = 0;

	for(i = 0; i < rec_count; ++i)
	{
		unsigned char *cur_quality = records[i].quality;
		for(j = 0; j < records[i].quality_len; ++j)
		{
			if(cur_quality[j] == prev && run_len < 254)
				++run_len;
			else
			{
				if(prev)
				{
					qua_stream[qua_idx++] = prev;
					run_stream[run_idx++] = (unsigned char) run_len;
				}

				run_len = 0;
				prev = cur_quality[j];
			}
		}
	}
	qua_stream[qua_idx++] = prev;
	run_stream[run_idx++] = (unsigned char) run_len;

	qua_stream_len = qua_idx;
	run_stream_len = run_idx;
	
	return true;
}

// ********************************************************************************************
bool CBlock::MakeUnRLE()
{
	uint32 i, j;
	uint32 tot_qua_len = 0;
	
	for(i = 0; i < rec_count; ++i)
		tot_qua_len += records[i].quality_len;

	uint32 run_len = 0;
	unsigned char qua = 0;
	uint32 qua_idx = 0;
	for(i = 0; i < rec_count; ++i)
	{
		no_of_N[i] = 0;
		unsigned char *cur_quality = records[i].quality;
		for(j = 0; j < records[i].quality_len; ++j)
		{
			if(run_len == 0)
			{
				qua = qua_stream[qua_idx];
				run_len = run_stream[qua_idx]+1;
				qua_idx++;
			}
			records[i].quality[j] = qua;
			--run_len;	
			if(qua >= 128)
				no_of_N[i]++;
		}
	}

	return true;
}

