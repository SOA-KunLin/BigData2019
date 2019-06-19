## Usage

Build docker image

```shell
docker build -t <hub-user>/<repo-name> .
docker push <hub-user>/<repo-name>
```

Prepare reference

```shell
wget ftp://ftp.ncbi.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/special_requests/GRCh37-lite.fa.gz
gunzip GRCh37-lite.fa.gz
mkdir final_project/data/Reference
mv GRCh37-lite.fa final_project/data/Reference/GRCh37-lite.fa
bwa index final_project/data/Reference/GRCh37-lite.fa
samtools faidx final_project/data/Reference/GRCh37-lite.fa
```

Run cromwell

```shell
java -jar -Dconfig.file=application.conf <path/to/cromwell>/cromwell-<version>.jar run bigdata.wdl --options options.json --inputs inputs.json
```
