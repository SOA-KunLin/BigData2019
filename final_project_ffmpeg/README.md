## Usage

Download Cromwell

[https://software.broadinstitute.org/wdl/](https://software.broadinstitute.org/wdl/)

Install Docker

[https://docs.docker.com/install/linux/docker-ce/ubuntu/](https://docs.docker.com/install/linux/docker-ce/ubuntu/)

Install Google Cloud SDK

[https://cloud.google.com/sdk/downloads#apt-get](https://cloud.google.com/sdk/downloads#apt-get)

Build Docker image

```shell
docker build -t <hub-user>/<repo-name> .
docker push <hub-user>/<repo-name>
```

Upload material files to Google Cloud Storage

```shell
gsutil mb gs://yourbuckets
gsutil -m cp path/to/Charizard/frame_*.jpg gs://yourbuckets/Charizard/
gsutil -m cp path/to/Charizard/water_marker.png gs://yourbuckets/Charizard/
gsutil -m cp path/to/Charizard/audio.mp3 gs://yourbuckets/Charizard/
```

Run Cromwell

```shell
java -Dconfig.file=google.conf -jar <path/to/cromwell>/cromwell-<version>.jar run ffmpeg.wdl --options options.json --inputs inputs.json
```
