#!/bin/bash
for i in {0..99}
do
    printf -v groupname "groups_%03d" $i
    mkdir $groupname
    cd $groupname
    for j in {0..447}
    do
        wget --content-disposition --header="API-Key: <enter API key here>" "http://www.tng-project.org/api/TNG100-1/files/groupcat-$i.$j.hdf5"
    done
    cd ..
done
