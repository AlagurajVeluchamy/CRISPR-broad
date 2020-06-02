#! /usr/bin/python
from argumentparse import *
from mapbwa import *

### Needed modules biopython

##############################################################################
#                               MAIN
##############################################################################
if __name__ == "__main__":
    arguments = arg_parsing()
    file = arguments.func(arguments)
    mapbwafastatogenome(arguments)
    filtersam(arguments)
    overlapeachchromosome(arguments)
    print file
