Code to replicate the HPRC and HGSVC experiments in "IsoformCheck: Protein isoform analysis from de novo genome assemblies."

To run the experiments you need to have python and R installed. Then:
- `git submodule update --init --recursive`
- `wget https://zenodo.org/records/17541941/files/IsoformCheck_db_v1_20251015.tar.gz?download=1 -O IsoformCheck_db_v1_20251015.tar.gz`
- `tar -xzf IsoformCheck_db_v1_20251015.tar.gz`
- `sh script.sh`
