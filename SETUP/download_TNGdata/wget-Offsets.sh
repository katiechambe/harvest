#!/bin/bash
for i in {0..99}
do
    wget --content-disposition --header="API-Key: <enter API key here>" "http://www.tng-project.org/api/TNG100-1/files/offsets.$i.hdf5"
done

