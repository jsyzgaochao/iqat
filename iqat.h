#ifndef _IQAT_H
#define _IQAT_H

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>

extern "C"
{
#include "libavutil/avutil.h"
#include "libavcodec/avcodec.h"
#include "libavformat/avformat.h"
#include "libswscale/swscale.h"
}

typedef struct
{
    char               filepath[256];
    bool               israw;
    FILE              *fp;
    AVFormatContext   *pFormatCtx;
    int                videoStream;
    AVCodecContext    *pCodecCtxOrig;
    AVCodecContext    *pCodecCtx;
    AVCodec           *pCodec;
    AVFrame           *pFrame;
    AVFrame           *pFrameScale;
    AVPacket           packet;
    bool               isscale;
    struct SwsContext *sws_ctx;
    uint8_t           *scale_buffer;
    uint8_t *y_data;
    uint8_t *u_data;
    uint8_t *v_data;
    double *y_psnr, *u_psnr, *v_psnr, *psnr, *ssim;
    double psnr_mean, y_psnr_mean, u_psnr_mean, v_psnr_mean, ssim_mean;
    double psnr_stdv, ssim_stdv;
} VideoInfo;

typedef struct
{
    int width;
    int height;
    bool log;
    VideoInfo *ref_video;
    VideoInfo *video[256];
    int video_num;
    char *last_error_info;
    int err_code;
    int frames;
} Params;

enum
{
    ERR_NONE = 0,
    ERR_HELP = -1,
    ERR_UNKNOWN_OPTION = -2,
    ERR_WRONG_FILE = -3,
    ERR_WRONG_PARAMS = -4,
    ERR_FFMPEG = -5,
    ERR_NO_MEMORY = -6
};

#define ERR_RET(err_code, params) { set_err_code(err_code, params); return err_code; }
#define x264_alloca(x) (void*)(((intptr_t)alloca((x)+15)+15)&~15)
#define XCHG(type,a,b) do { type t = a; a = b; b = t; } while( 0 )
#define X264_MIN(a,b) ( (a)<(b) ? (a) : (b) )
#ifdef _WIN32
#define strcasecmp stricmp
#endif

void print_help(char *strAppName, const char *strErrorMessage, ...);
int check_inputfile(Params *params);
void set_err_code(int errcode, Params* params);
int check_params(Params* params);
VideoInfo *new_videoinfo();
int parse(char* argv[], int argc, Params* params);
int init_video(VideoInfo *vi, Params* params);
int get_frame(VideoInfo *vi, Params* params);
void free_video(VideoInfo *vi);

void ssim_4x4x2_core(const uint8_t *pix1, int stride1,
                     const uint8_t *pix2, int stride2,
                     int sums[2][4]);
float ssim_end1(int s1, int s2, int ss, int s12);
float ssim_end4(int sum0[5][4], int sum1[5][4], int width);
float x264_pixel_ssim_wxh(
        uint8_t *pix1, int stride1,
        uint8_t *pix2, int stride2,
        int width, int height);

#endif

