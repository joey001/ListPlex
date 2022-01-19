# ListPlex
A fast k-plex listing algorithm **Without Prescribed Size (Pro1)**

# Setup
```shell
make
```

# Usage
  ./listPlex \<file\> \<k\> [thread-number]

e.g.

```bash
./listPlex ./jazz.bin 3 1
```

# Format
The input graph should be a binary format.
One can convert an edge list format graph file (SNAP format) into this binary format by a converter `toBin` .

usage:

  ./toBin \<input\> \<output\>

e.g.

```bash
./toBin ./jazz.txt ./jazz.bin
```

