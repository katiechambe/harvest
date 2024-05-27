File structure for working with data:

```
IllustrisTNG/ (parent directory)
├── output/
│   ├── groups_***/
│   |    ├── fof_subhalo_tab_0**.***.hdf5
|
├── postprocessing/
│   ├── offsets/
│   |    ├── offsets_***.hdf5
|   ├── tree_extended.**.hdf5
|
├── trees/



- output/ (contains group catalogs)
    - groups_000
    - ...
    - groups_099
 
- 
    - offsets/
        - offsets_000.hdf5
        - ...
        - offsets_099.hdf5

- trees/
```


https://www.tng-project.org/data/docs/scripts/