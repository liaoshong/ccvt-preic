## Generating pre-initial conditions for cosmological *N*-body simulations with capacity constrained Voronoi tessellations

Copyright (C) 2018  Shihong Liao

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

---

Following is a quick summary to get started with using the code. Please
refer to **User Guide** (*doc/user-guide.pdf*) for further details.

## Downloads
You can directly download this code [here](https://github.com/liaoshong/ccvt-preic),
or use git to clone it:
```
$ git clone https://github.com/liaoshong/ccvt-preic.git
```

## Environment requirements
To compile the program, you need an ANSI C compiler with OpenMP
supported (e.g. [GCC](https://gcc.gnu.org/) later than version 4.2).

## Compile
To compile the code, type
```
$ make
```
in the command line, and normally it will generate an executable named
*ccvt-preic*.

## Run
Run the code with the below command.
```
$ ./ccvt-preic results/run.param
```

## Outputs
The code outputs two files:

1. *ccvt_particle_[NUM_PART_EACH_DIM]_capacity_[NUM_CAPACITY_EACH_DIM].txt*: a
text file with particles' positions.
2. *gadget_particle_[NUM_PART_EACH_DIM]_capacity_[NUM_CAPACITY_EACH_DIM]*: a
file with the GADGET format (i.e., SnapFormat = 1 in GADGET's conventions).
This file can be directly used by the N-GenIC or 2LPTic codes (see User
Guide for details).

## Contacts
Questions and bug reports can be sent to [liaoshong@gmail.com](mailto:liaoshong@gmail.com).
Feedbacks and comments are welcome.
