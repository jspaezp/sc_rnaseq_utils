#!/bin/bash

# Takes a .dat file and makes a lot of plots rotating it ...
set -e

if [[ -f $@ ]] 
then
  echo "Using ${@} as data file"
else 
  echo "Please pass a file as an argument"
fi

CURRDIR=$(pwd)
TEMPDIR=$(mktemp -d)
DATAPATH=$(readlink -f "${@}")
trap "rm -rf ${TEMPDIR}" EXIT

echo "Changing Directory to: ${TEMPDIR}"
cd "${TEMPDIR}"

echo "Writting palette file"
cat > palette.pal << EOF
# line styles for ColorBrewer Spectral
# for use with divering data
# provides 8 rainbow colors with red low, pale yellow middle, and blue high
# compatible with gnuplot >=4.2
# author: Anna Schneider

# line styles
set style line 1 lt 1 lc rgb '#D53E4F' # red
set style line 2 lt 1 lc rgb '#F46D43' # orange
set style line 3 lt 1 lc rgb '#FDAE61' # pale orange
set style line 4 lt 1 lc rgb '#FEE08B' # pale yellow-orange
set style line 5 lt 1 lc rgb '#E6F598' # pale yellow-green
set style line 6 lt 1 lc rgb '#ABDDA4' # pale green
set style line 7 lt 1 lc rgb '#66C2A5' # green
set style line 8 lt 1 lc rgb '#3288BD' # blue

# palette
set palette defined ( 0 '#D53E4F',\
          1 '#F46D43',\
		      2 '#FDAE61',\
		      3 '#FEE08B',\
		      4 '#E6F598',\
		      5 '#ABDDA4',\
		      6 '#66C2A5',\
		      7 '#3288BD' )
EOF

echo "Simlinking data file"
ln -s "${DATAPATH}" . --verbose

pwd

echo "Plotting"
FULLCOMMAND="load 'palette.pal' ;\
set datafile separator \",\" ;\
unset xtics ;\
unset ytics ;\
unset ztics ;\
unset border  ;\
set pointsize 0.5 ;\
do for [rot=1:59] { ;\
    set term png        ;\
    outfile = sprintf('${TEMPDIR}/foo%03.0f.png',rot) ;\
    set output outfile ;\
    set view 60,6*rot,1,1 ;\
    splot \"${DATAPATH}\" with points pointtype 7 palette ;\
    replot ;\
    set term x11  ;\
}"
gnuplot -p -e "${FULLCOMMAND}" 

echo "Found $(ls -1 *.png | wc -l) PNG Files"

echo "Writting final GIF as: ${DATAPATH}.gif"
convert -loop 0 *.png "${DATAPATH}.gif"

echo "Cleaning temporary Directory"
rm -rf ${TEMPDIR}
