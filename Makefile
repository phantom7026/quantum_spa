
TRGT=SPA

OBJS=main_rp.o channel.o utility.o syndrome.o galois_field.o
SRCS=$(OBJS:.o=.cpp)

$(TRGT): $(OBJS)
	mpicxx -Ofast -m64 -o $@ $^ 
	
.cpp.o:
	mpic++ -Ofast -g -c $<

deep:
	gccmakedep $(SRCS)
	
clean:
	rm -rf $(OBJS) $(TRGT) *.o
	
