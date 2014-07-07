#! /bin/sh

if [ $# -lt 1 ]; then
  echo 1>&2 Usage: $0 "(runA | runB | runC)" [number_of_threads]
  exit -1
fi

if [ ! -f ../$1/input/BarnesHut.in ]; then
  echo 1>&2 Error: cannot find input file ../$1/input/BarnesHut.in
  echo 1>&2 Usage: $0 "(runA | runB | runC)" number_of_threads
  exit -1
fi

if [ ! -f ./BarnesHut ]; then
  DIR=`pwd`
  cd .. && make all && make install
  cd $DIR
fi

./BarnesHut ../$1/input/BarnesHut.in $2 > /dev/null
