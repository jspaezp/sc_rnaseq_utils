#!/bin/bash

SCRIPTLOCATION=$(readlink -f ./gp_3dscat.bash)
TEMPDIR=$(mktemp -d)
trap "rm -rf ${TEMPDIR}" EXIT
cd "${TEMPDIR}"

cat > data.dat << EOF
#x,y,z,col
1,2,3,4
2,2,3,4
2,3,4,5
3,3,4,5
EOF

bash "${SCRIPTLOCATION}" data.dat

if [[ -f data.dat.gif ]]
then
  echo "Successfully generates a gif"
else
  echo "Does not generate gif"
  EXIT 1
fi

rm -rf ${TEMPDIR}

