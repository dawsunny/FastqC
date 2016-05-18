all: dsrc

CFLAGS  = -O3 -m64 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE

.cpp.o:
	g++ $(CFLAGS) -c $< 

dsrc: dsrc.o block.o compress.o DSRCFile.o field.o huffman.o io.o lz.o superblock.o
	g++ $(CFLAGS) -o $@ dsrc.o block.o compress.o DSRCFile.o field.o huffman.o io.o lz.o superblock.o

clean:
	-rm *.o
	-rm dsrc
