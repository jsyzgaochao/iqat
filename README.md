iqat (Image Quality Assessment Tool)
==========
It is a useful tool for video image quality assessment, which can calculate the PSNR and SSIM of multi-channel videos compared to the reference video. The video format can be either raw data or compressed format (such as h.265, h.264, mpeg4, etc).

Requirements:
--------
- FFmpeg
- OpenMP
- SSE, AVX instruction set

FFmpeg helps to decode the compressed format and resize the image. OpenMP helps to execute parallelly for multi videos. The processor should support SSE, AVX instruction set for SIMD vectorization.

Running:
--------
```
Usage: iqat -r RefFile File1 [File2 ...] [options]

Options:
   [-w width]         - video width
   [-h height]        - video height
   [-n num]           - frames number
   [-skip num]        - frame interval
   [-log]             - output log files

Example: iqat -r RefVideo.xxx Video1.xxx Video2.xxx -log -w 1280 -h 720 -n 1000 -skip 2
```
If RefVideo.xxx is a compressed video containing resolution information, the width and height can be omitted in the command line. -log option will generate a log file for each video, which contains image quality information of each frame. -n option specifies the numbers of frame to test. -skip option specifies the interval of frames to test.
