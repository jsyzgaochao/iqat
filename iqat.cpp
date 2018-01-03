#include <omp.h>
#include "iqat.h"

unsigned long get_msec()
{
    struct timeval tval;
    gettimeofday(&tval, NULL);
    return tval.tv_sec * 1000 + (tval.tv_usec + 500) / 1000;
}

void print_help(char *strAppName, const char *strErrorMessage, ...)
{
    if (strErrorMessage)
    {
        va_list args;
        printf("ERROR: ");
        va_start(args, strErrorMessage);
        vprintf(strErrorMessage, args);
        va_end(args);
        printf("\n\n");
    }
    printf("Image Quality Assessment Tool by Gao Chao\n");
    printf("Usage: %s -r RefFile File1 [File2 ...] [options]\n", strAppName);
    printf("\n");
    printf("Options:\n");
    printf("   [-w width]         - video width\n");
    printf("   [-h height]        - video height\n");
    printf("   [-n num]           - frame numbers\n");
    printf("   [-log]             - output log files\n");
}

int check_inputfile(Params *params)
{
    FILE *fp = fopen(params->ref_video->filepath, "r");
    if (NULL == fp)
    {
        params->last_error_info = params->ref_video->filepath;
        return ERR_WRONG_FILE;
    }
    else
        fclose(fp);
    for (int i = 0; i < params->video_num; i++)
    {
        FILE *fp = fopen(params->video[i]->filepath, "r");
        if (NULL == fp)
        {
            params->last_error_info = params->video[i]->filepath;
            return ERR_WRONG_FILE;
        }
        else
            fclose(fp);
    }
    return ERR_NONE;
}

void set_err_code(int errcode, Params* params)
{
    if(ERR_NONE == params->err_code)
        params->err_code = errcode;
}

int check_params(Params* params)
{
    if (NULL == params->ref_video)
    {
        params->last_error_info = "No Reference Video !";
        return ERR_WRONG_PARAMS;
    }
    if (0 == params->video_num)
    {
        params->last_error_info = "No Videos !";
        return ERR_WRONG_PARAMS;
    }

    params->ref_video->israw = false;
    char *ext = strrchr(params->ref_video->filepath, '.');
    if (ext)
        ext++;

    if ((0 == strcasecmp(ext, "yuv")) || (0 == strcasecmp(ext, "raw")))
        params->ref_video->israw = true;

    if (params->ref_video->israw)
    {
        if (params->width == 0 || params->height == 0)
        {
            params->last_error_info = "Invalid Width and Height for Raw Video !";
            return ERR_WRONG_PARAMS;
        }
    }
    else
    {
        if (params->width == 0 || params->height == 0)
        {
            AVFormatContext   *pFormatCtx = NULL;
            int                videoStream;
            AVCodecContext    *pCodecCtx = NULL;

            // Open video file
            if (avformat_open_input(&pFormatCtx, params->ref_video->filepath, NULL, NULL) != 0)
                return ERR_FFMPEG; // Couldn't open file

            // Retrieve stream information
            if (avformat_find_stream_info(pFormatCtx, NULL)<0)
                return ERR_FFMPEG; // Couldn't find stream information

            // Find the first video stream
            videoStream = -1;
            for (int i = 0; i<pFormatCtx->nb_streams; i++)
                if (pFormatCtx->streams[i]->codec->codec_type == AVMEDIA_TYPE_VIDEO) {
                    videoStream = i;
                    break;
                }
            if (videoStream == -1)
                return ERR_FFMPEG; // Didn't find a video stream

            // Get a pointer to the codec context for the video stream
            pCodecCtx = pFormatCtx->streams[videoStream]->codec;

            params->width = pCodecCtx->width;
            params->height = pCodecCtx->height;

            // Close the codecs
            avcodec_close(pCodecCtx);

            // Close the video file
            avformat_close_input(&pFormatCtx);
        }
    }

    return ERR_NONE;
}

VideoInfo *new_videoinfo()
{
    VideoInfo *vi = (VideoInfo*)malloc(sizeof(VideoInfo));
    memset(vi, 0, sizeof(VideoInfo));
    return vi;
}

int parse(char* argv[], int argc, Params* params)
{
    if (1 == argc)
    {
        print_help(argv[0], NULL);
        ERR_RET(ERR_HELP,params);
    }
    params->video_num = 0;
    params->log = false;
    for (int i = 1; i < argc; i++)
    {
        if ('-' != argv[i][0])
        {
            params->video[params->video_num] = new_videoinfo();
            strcpy(params->video[params->video_num]->filepath, argv[i]);
            params->video_num++;
        }
        else if (0 == strcmp(argv[i], "-?"))
        {
            print_help(argv[0], NULL);
            ERR_RET(ERR_HELP, params);
        }
        else if (0 == strcmp(argv[i], "-w"))
        {
            params->width = atoi(argv[++i]);
        }
        else if (0 == strcmp(argv[i], "-h"))
        {
            params->height = atoi(argv[++i]);
        }
        else if (0 == strcmp(argv[i], "-r"))
        {
            params->ref_video = new_videoinfo();
            strcpy(params->ref_video->filepath, argv[++i]);
        }
        else if (0 == strcmp(argv[i], "-n"))
        {
            params->frames_num = atoi(argv[++i]);
        }
        else if (0 == strcmp(argv[i], "-log"))
        {
            params->log = true;
        }
        else
        {
            print_help(argv[0], "Unknown options: %s", argv[i]);
            ERR_RET(ERR_UNKNOWN_OPTION, params);
        }
    }
    if (check_inputfile(params))
    {
        print_help(argv[0], "Wrong input file : %s", params->last_error_info);
        ERR_RET(ERR_WRONG_FILE, params);
    }
    if (check_params(params))
    {
        print_help(argv[0], "Wrong params : %s", params->last_error_info);
        ERR_RET(ERR_WRONG_PARAMS, params);
    }

    return ERR_NONE;
}

int init_video(VideoInfo *vi, Params* params)
{
    vi->israw = false;
    char *ext = strrchr(vi->filepath, '.');
    if (ext)
        ext++;

    if ((0 == strcasecmp(ext, "yuv")) || (0 == strcasecmp(ext, "raw")))
        params->ref_video->israw = true;

    if(vi->israw)
    {
        vi->fp = fopen(vi->filepath, "rb");
    }
    else
    {
        // Open video file
        if (avformat_open_input(&vi->pFormatCtx, vi->filepath, NULL, NULL) != 0)
            return ERR_FFMPEG; // Couldn't open file

        // Retrieve stream information
        if (avformat_find_stream_info(vi->pFormatCtx, NULL)<0)
            return ERR_FFMPEG; // Couldn't find stream information

        // Dump information about file onto standard error
        av_dump_format(vi->pFormatCtx, 0, vi->filepath, 0);

        // Find the first video stream
        vi->videoStream = -1;
        for (int i = 0; i<vi->pFormatCtx->nb_streams; i++)
            if (vi->pFormatCtx->streams[i]->codec->codec_type == AVMEDIA_TYPE_VIDEO) {
                vi->videoStream = i;
                break;
            }
        if (vi->videoStream == -1)
            return ERR_FFMPEG; // Didn't find a video stream

        // Get a pointer to the codec context for the video stream
        vi->pCodecCtxOrig = vi->pFormatCtx->streams[vi->videoStream]->codec;
        // Find the decoder for the video stream
        vi->pCodec = avcodec_find_decoder(vi->pCodecCtxOrig->codec_id);
        if (vi->pCodec == NULL) {
            fprintf(stderr, "Unsupported codec!\n");
            return ERR_FFMPEG; // Codec not found
        }
        // Copy context
        vi->pCodecCtx = avcodec_alloc_context3(vi->pCodec);
        if (avcodec_copy_context(vi->pCodecCtx, vi->pCodecCtxOrig) != 0) {
            fprintf(stderr, "Couldn't copy codec context");
            return ERR_FFMPEG; // Error copying codec context
        }

        // Open codec
        if (avcodec_open2(vi->pCodecCtx, vi->pCodec, NULL)<0)
            return ERR_FFMPEG; // Could not open codec

        // Allocate video frame
        vi->pFrame = av_frame_alloc();

        if (vi->pCodecCtx->width == params->width && vi->pCodecCtx->height == params->height)
            vi->isscale = false;
        else
            vi->isscale = true;

        if (vi->isscale)
        {
            // Allocate an AVFrame structure
            vi->pFrameScale = av_frame_alloc();
            if (vi->pFrameScale == NULL)
                return ERR_FFMPEG;

            // Determine required buffer size and allocate buffer
            int numBytes = avpicture_get_size(AV_PIX_FMT_YUV420P, params->width,
                                              params->height);
            vi->scale_buffer = (uint8_t *)av_malloc(numBytes * sizeof(uint8_t));

            // Assign appropriate parts of buffer to image planes in pFrameRGB
            // Note that pFrameRGB is an AVFrame, but AVFrame is a superset
            // of AVPicture
            avpicture_fill((AVPicture *)vi->pFrameScale, vi->scale_buffer, AV_PIX_FMT_YUV420P,
                           params->width, params->height);

            // initialize SWS context for software scaling
            vi->sws_ctx = sws_getContext(vi->pCodecCtx->width,
                                         vi->pCodecCtx->height,
                                         vi->pCodecCtx->pix_fmt,
                                         params->width,
                                         params->height,
                                         AV_PIX_FMT_YUV420P,
                                         SWS_BICUBIC,
                                         NULL,
                                         NULL,
                                         NULL
                                         );
        }
    }

    vi->y_data = (uint8_t *)malloc(sizeof(uint8_t)*params->width*params->height);
    vi->u_data = (uint8_t *)malloc(sizeof(uint8_t)*params->width*params->height / 4);
    vi->v_data = (uint8_t *)malloc(sizeof(uint8_t)*params->width*params->height / 4);
    vi->size = 0;

    return ERR_NONE;
}

int get_frame(VideoInfo *vi, Params* params)
{
    if (vi->israw)
    {
        if (1 != fread(vi->y_data, sizeof(uint8_t)*params->width*params->height, 1, vi->fp)) return 0;
        if (1 != fread(vi->u_data, sizeof(uint8_t)*params->width*params->height / 4, 1, vi->fp)) return 0;
        if (1 != fread(vi->v_data, sizeof(uint8_t)*params->width*params->height / 4, 1, vi->fp)) return 0;
        vi->size += params->width*params->height*3/2;
        return 1;
    }
    else
    {
        int frameFinished;
        while (av_read_frame(vi->pFormatCtx, &vi->packet) >= 0) {
            // Is this a packet from the video stream?
            if (vi->packet.stream_index == vi->videoStream) {
                // Decode video frame
                avcodec_decode_video2(vi->pCodecCtx, vi->pFrame, &frameFinished, &vi->packet);

                // Did we get a video frame?
                if (frameFinished) {
                    if (vi->isscale)
                    {
                        sws_scale(vi->sws_ctx, vi->pFrame->data,
                                  vi->pFrame->linesize, 0, vi->pCodecCtx->height,
                                  vi->pFrameScale->data, vi->pFrameScale->linesize);
                        memcpy(vi->y_data, vi->pFrameScale->data[0], sizeof(uint8_t)*params->width*params->height);
                        memcpy(vi->u_data, vi->pFrameScale->data[1], sizeof(uint8_t)*params->width*params->height / 4);
                        memcpy(vi->v_data, vi->pFrameScale->data[2], sizeof(uint8_t)*params->width*params->height / 4);
                    }
                    else
                    {
                        memcpy(vi->y_data, vi->pFrame->data[0], sizeof(uint8_t)*params->width*params->height);
                        memcpy(vi->u_data, vi->pFrame->data[1], sizeof(uint8_t)*params->width*params->height / 4);
                        memcpy(vi->v_data, vi->pFrame->data[2], sizeof(uint8_t)*params->width*params->height / 4);
                    }
                    vi->size += vi->packet.size;
                    av_free_packet(&vi->packet);
                    return 1;
                }
            }
            av_free_packet(&vi->packet);
        }
        return 0;
    }
}

void free_video(VideoInfo *vi)
{
    free(vi->y_data);
    free(vi->u_data);
    free(vi->v_data);
    free(vi->psnr);
    free(vi->y_psnr);
    free(vi->u_psnr);
    free(vi->v_psnr);
    free(vi->ssim);

    if (vi->israw)
    {
        fclose(vi->fp);
    }
    else
    {
        if (vi->isscale)
        {
            // Free the scaled image
            sws_freeContext(vi->sws_ctx);
            av_free(vi->scale_buffer);
            av_frame_free(&vi->pFrameScale);
        }

        // Free the YUV frame
        av_frame_free(&vi->pFrame);

        // Close the codec
        avcodec_close(vi->pCodecCtx);
        avcodec_close(vi->pCodecCtxOrig);

        // Close the video file
        avformat_close_input(&vi->pFormatCtx);
    }
}

void ssim_4x4x2_core(const uint8_t *pix1, int stride1,
                     const uint8_t *pix2, int stride2,
                     int sums[2][4])
{
    int x, y, z;
    for (z = 0; z<2; z++)
    {
        uint32_t s1 = 0, s2 = 0, ss = 0, s12 = 0;
        for (y = 0; y<4; y++)
            for (x = 0; x<4; x++)
            {
                int a = pix1[x + y*stride1];
                int b = pix2[x + y*stride2];
                s1 += a;
                s2 += b;
                ss += a*a;
                ss += b*b;
                s12 += a*b;
            }
        sums[z][0] = s1;
        sums[z][1] = s2;
        sums[z][2] = ss;
        sums[z][3] = s12;
        pix1 += 4;
        pix2 += 4;
    }
}

float ssim_end1(int s1, int s2, int ss, int s12)
{
    static const int ssim_c1 = (int)(.01*.01 * 255 * 255 * 64 + .5);
    static const int ssim_c2 = (int)(.03*.03 * 255 * 255 * 64 * 63 + .5);
    int vars = ss * 64 - s1*s1 - s2*s2;
    int covar = s12 * 64 - s1*s2;
    return (float)(2 * s1*s2 + ssim_c1) * (float)(2 * covar + ssim_c2)\
            / ((float)(s1*s1 + s2*s2 + ssim_c1) * (float)(vars + ssim_c2));
}

float ssim_end4(int sum0[5][4], int sum1[5][4], int width)
{
    int i;
    float ssim = 0.0;
    for (i = 0; i < width; i++)
        ssim += ssim_end1(sum0[i][0] + sum0[i + 1][0] + sum1[i][0] + sum1[i + 1][0],
                sum0[i][1] + sum0[i + 1][1] + sum1[i][1] + sum1[i + 1][1],
                sum0[i][2] + sum0[i + 1][2] + sum1[i][2] + sum1[i + 1][2],
                sum0[i][3] + sum0[i + 1][3] + sum1[i][3] + sum1[i + 1][3]);
    return ssim;
}

float x264_pixel_ssim_wxh(
        uint8_t *pix1, int stride1,
        uint8_t *pix2, int stride2,
        int width, int height)
{
    int x, y, z;
    float ssim = 0.0;
    int(*sum0)[4] = (int(*)[4])x264_alloca(4 * (width / 4 + 3) * sizeof(int));
    int(*sum1)[4] = (int(*)[4])x264_alloca(4 * (width / 4 + 3) * sizeof(int));
    width >>= 2;
    height >>= 2;
    z = 0;
    for (y = 1; y < height; y++)
    {
        for (; z <= y; z++)
        {
            XCHG( auto, sum0, sum1 );
            for (x = 0; x < width; x += 2)
                ssim_4x4x2_core(&pix1[4 * (x + z*stride1)], stride1, &pix2[4 * (x + z*stride2)], stride2, &sum0[x]);
        }
        for (x = 0; x < width - 1; x += 4)
            ssim += ssim_end4(sum0 + x, sum1 + x, X264_MIN(4, width - x - 1));
    }
    return ssim / ((width - 1) * (height - 1));
}

int main(int argc, char *argv[])
{
    av_register_all();
    Params *params = (Params*)malloc(sizeof(Params));
    memset(params, 0, sizeof(Params));
    parse(argv, argc, params);
    if (0 != params->err_code)
        return params->err_code;

    init_video(params->ref_video, params);
    for (int i = 0; i < params->video_num; i++)
    {
        init_video(params->video[i], params);
    }

    int frames = 0;
    int buffer_size = 0;
    int y_size = params->width*params->height;
    int uv_size = params->width*params->height / 4;
    int total_size = y_size + 2 * uv_size;
    long last_msec = get_msec();

    printf("\n");

    while (get_frame(params->ref_video, params))
    {
        int break_flag = 0;
        if (frames >= buffer_size)
        {
            buffer_size += 0xffff;
            for (int i = 0; i < params->video_num; i++)
            {
                if (!(params->video[i]->psnr = (double*)realloc(params->video[i]->psnr, buffer_size * sizeof(double))))
                {
                    printf("ERROR: Not enough memory\n\n");
                    return ERR_NO_MEMORY;
                }
                if (!(params->video[i]->y_psnr = (double*)realloc(params->video[i]->y_psnr, buffer_size * sizeof(double))))
                {
                    printf("ERROR: Not enough memory\n\n");
                    return ERR_NO_MEMORY;
                }
                if (!(params->video[i]->u_psnr = (double*)realloc(params->video[i]->u_psnr, buffer_size * sizeof(double))))
                {
                    printf("ERROR: Not enough memory\n\n");
                    return ERR_NO_MEMORY;
                }
                if (!(params->video[i]->v_psnr = (double*)realloc(params->video[i]->v_psnr, buffer_size * sizeof(double))))
                {
                    printf("ERROR: Not enough memory\n\n");
                    return ERR_NO_MEMORY;
                }
                if (!(params->video[i]->ssim = (double*)realloc(params->video[i]->ssim, buffer_size * sizeof(double))))
                {
                    printf("ERROR: Not enough memory\n\n");
                    return ERR_NO_MEMORY;
                }
            }
        }
#pragma omp parallel for
        for (int i = 0; i < params->video_num; i++)
        {
            if (get_frame(params->video[i], params))
            {
                double y_mse = 0.0, u_mse = 0.0, v_mse = 0.0, total_mse = 0.0;
                double diff = 0.0;
                for (int j = 0; j < y_size; j++)
                {
                    diff = params->video[i]->y_data[j] - params->ref_video->y_data[j];
                    y_mse += diff * diff;
                }
                for (int j = 0; j < uv_size; j++)
                {
                    diff = params->video[i]->u_data[j] - params->ref_video->u_data[j];
                    u_mse += diff * diff;
                    diff = params->video[i]->v_data[j] - params->ref_video->v_data[j];
                    v_mse += diff * diff;
                }
                total_mse = y_mse + u_mse + v_mse;
                params->video[i]->psnr[frames] = total_mse > 0.0 ? 20 * (log10(255 / sqrt(total_mse / total_size))) : 100;
                params->video[i]->y_psnr[frames] = y_mse > 0.0 ? 20 * (log10(255 / sqrt(y_mse / y_size))) : 100;
                params->video[i]->u_psnr[frames] = u_mse > 0.0 ? 20 * (log10(255 / sqrt(u_mse / uv_size))) : 100;
                params->video[i]->v_psnr[frames] = v_mse > 0.0 ? 20 * (log10(255 / sqrt(v_mse / uv_size))) : 100;
                params->video[i]->ssim[frames] = x264_pixel_ssim_wxh(params->ref_video->y_data, params->width, params->video[i]->y_data, params->width, params->width, params->height);
                params->video[i]->psnr_mean += params->video[i]->psnr[frames];
                params->video[i]->y_psnr_mean += params->video[i]->y_psnr[frames];
                params->video[i]->u_psnr_mean += params->video[i]->u_psnr[frames];
                params->video[i]->v_psnr_mean += params->video[i]->v_psnr[frames];
                params->video[i]->ssim_mean += params->video[i]->ssim[frames];
            }
            else
            {
                break_flag = 1;
            }
        }
        if (break_flag)
            break;
        frames++;
        long now_msec = get_msec();
        if (now_msec - last_msec > 1000)
        {
            printf("\rframes: %d", frames);
            fflush(stdout);
            last_msec = now_msec;
        }
        if (params->frames_num && frames == params->frames_num)
            break;
    }

    printf("\rtotal frames: %d\n\n", frames);
    fflush(stdout);

    params->frames = frames;
    if (params->frames > 0)
    {
#pragma omp parallel for
        for (int i = 0; i < params->video_num; i++)
        {
            params->video[i]->psnr_mean /= params->frames;
            params->video[i]->y_psnr_mean /= params->frames;
            params->video[i]->u_psnr_mean /= params->frames;
            params->video[i]->v_psnr_mean /= params->frames;
            params->video[i]->ssim_mean /= params->frames;

            double diff;
            for (int j = 0; j < params->frames; j++)
            {
                diff = params->video[i]->psnr[j] - params->video[i]->psnr_mean;
                params->video[i]->psnr_stdv += diff * diff;
                diff = params->video[i]->ssim[j] - params->video[i]->ssim_mean;
                params->video[i]->ssim_stdv += diff * diff;
            }
            params->video[i]->psnr_stdv = sqrt(params->video[i]->psnr_stdv / params->frames);
            params->video[i]->ssim_stdv = sqrt(params->video[i]->ssim_stdv / params->frames);
        }

    }

    for (int i = 0; i < params->video_num; i++)
    {
        printf("Video: %s\n", params->video[i]->filepath);
        printf("Size: %0.2f MB\n", (float)params->video[i]->size/1024.0f/1024.0f);
        printf("Frames: %d\n", params->frames);
        printf("PSNR: mean %.3f Y %.3f U %.3f V %.3f stdv %.3f\n", params->video[i]->psnr_mean,
               params->video[i]->y_psnr_mean, params->video[i]->u_psnr_mean, params->video[i]->v_psnr_mean, params->video[i]->psnr_stdv);
        printf("SSIM: mean %.4f (%.3fdB) stdv %.4f\n\n", params->video[i]->ssim_mean,
               params->video[i]->ssim_mean < 1.0 ? -10 * (log10(1 - params->video[i]->ssim_mean)) : 100,
               params->video[i]->ssim_stdv);
        if (params->log)
        {
            char log_filename[256];
            strcpy(log_filename, params->video[i]->filepath);
            strcat(log_filename, ".log");
            FILE *fp = fopen(log_filename, "w");
            if (NULL != fp)
            {
                fprintf(fp, "Video: %s\n", params->video[i]->filepath);
                fprintf(fp, "Size: %0.2f MB\n", (float)params->video[i]->size/1024.0f/1024.0f);
                fprintf(fp, "Frames: %d\n", params->frames);
                fprintf(fp, "PSNR: mean %.3f Y %.3f U %.3f V %.3f stdv %.3f\n", params->video[i]->psnr_mean,
                        params->video[i]->y_psnr_mean, params->video[i]->u_psnr_mean, params->video[i]->v_psnr_mean, params->video[i]->psnr_stdv);
                fprintf(fp, "SSIM: mean %.4f (%.3fdB) stdv %.4f\n\n", params->video[i]->ssim_mean,
                        params->video[i]->ssim_mean < 1.0 ? -10 * (log10(1 - params->video[i]->ssim_mean)) : 100,
                        params->video[i]->ssim_stdv);
                fprintf(fp, "Frame\tPSNR\tY_PSNR\tU_PSNR\tV_PSNR\tSSIM\n");
                for (int j = 0; j < params->frames; j++)
                {
                    fprintf(fp, "%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", j + 1, params->video[i]->psnr[j], params->video[i]->y_psnr[j], params->video[i]->u_psnr[j],
                            params->video[i]->v_psnr[j], params->video[i]->ssim[j]);
                }
                fclose(fp);
            }
        }
    }

    free_video(params->ref_video);
    free(params->ref_video);
    for (int i = 0; i < params->video_num; i++)
    {
        free_video(params->video[i]);
        free(params->video[i]);
    }

    free(params);
}
