workflow myWorkflow {
    Array[File] images
    File water_marker
    File audio
    scatter (image in images) {
        call TaskAddWaterMarker { input:
	    image=image,
            water_marker=water_marker,
        }
    }
    call TaskMakeGIF { input: images=TaskAddWaterMarker.jpg }
    call TaskMakeMp4 { input: images=TaskAddWaterMarker.jpg }
    call TaskMakeAVI { input: images=TaskAddWaterMarker.jpg }
    call TaskMergeAudio { input:
        video=TaskMakeMp4.mp4,
        audio=audio,
    }
    output {
        File gif = TaskMakeGIF.gif
        File mp4 = TaskMakeMp4.mp4
        File AVI = TaskMakeAVI.avi
        File merge = TaskMergeAudio.mp4
    }
}

task TaskAddWaterMarker {
    File image
    File water_marker
    command {
        ffmpeg -y -i ${image} -i ${water_marker} -filter_complex "overlay=main_w-overlay_w:0" out.jpg
    }
    output {
        File jpg = "out.jpg"
    }
}

task TaskMakeGIF {
    Array[File] images
    command {
        i=1
        for fn in ${sep=" " images}; do
            cp $fn $(printf "frame_%02d.jpg" "$i")
            ((i++))
        done
        ffmpeg -y -i frame_%2d.jpg video.gif
    }
    output {
        File gif = "video.gif"
    }
}

task TaskMakeMp4 {
    Array[File] images
    command {
        i=1
        for fn in ${sep=" " images}; do
            cp $fn $(printf "frame_%02d.jpg" "$i")
            ((i++))
        done
        ffmpeg -y -i frame_%2d.jpg video.mp4
    }
    output {
        File mp4 = "video.mp4"
    }
}

task TaskMakeAVI {
    Array[File] images
    command {
        i=1
        for fn in ${sep=" " images}; do
            cp $fn $(printf "frame_%02d.jpg" "$i")
            ((i++))
        done
        ffmpeg -y -i frame_%2d.jpg video.avi
    }
    output {
        File avi = "video.avi"
    }
}

task TaskMergeAudio {
    File video
    File audio
    command {
        ffmpeg -y -i ${video} -i ${audio} -c:v copy -c:a aac -strict experimental -shortest merge.mp4
    }
    output {
        File mp4 = "merge.mp4"
    }
}
