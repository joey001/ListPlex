# ListPlex
A fast k-plex listing algorithm **With Prescribed Size (Pro2)**

# Setup
```shell
make
```

# Usage
  ./listPlex \<file\> \<k\> \<lb\> [thread-number]

e.g.

```bash
./listPlex ./jazz.bin 4 12 1
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

