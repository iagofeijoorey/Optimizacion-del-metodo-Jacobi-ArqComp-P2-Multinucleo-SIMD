
all:
	$(MAKE) -C src all

run:
	$(MAKE) -C src run ARGS="$(ARGS)"

clean:
	$(MAKE) -C src clean