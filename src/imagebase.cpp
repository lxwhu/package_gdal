// imagebase.cpp : Defines the exported functions for the DLL application.
//

#include "LxString.hpp"
#include "imagebase.h"
//#include <vector>
#include <gdal_priv.h>
#include "wavelet.h"

#ifndef WIN32
#define stricmp strcasecmp
#endif

struct image_dll_register
{
    image_dll_register()
    {
        GDALAllRegister();
        CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
    }
} g_image_dll_register;

#define  GEOTRANSFROM_CORNER2CENTER(adfGeoTransform) ApplyGeoTransform(adfGeoTransform, 0.5, 0.5, adfGeoTransform + 0, adfGeoTransform + 3)
#define  GEOTRANSFROM_CENTER2CORNER(adfGeoTransform) ApplyGeoTransform(adfGeoTransform, -0.5, -0.5, adfGeoTransform + 0, adfGeoTransform + 3)

#define BGR2GRAY(b,g,r)  ((r*38 + g*75 + b*15) >> 7)

class image
{
public:
    image(){    
        m_pImgSet = NULL;
        m_pMemDriver = NULL;
        m_imgCoor = left_top;
        m_tab16to8 = NULL;
        m_strImagePath[0] = 0;
    }
    ~image(){
        Reset();
        if (m_pImgSet) delete m_pImgSet;    m_pImgSet = NULL;
        if (m_tab16to8) delete[] m_tab16to8;    m_tab16to8 = NULL;
    }
    enum SAMP_STOR{ SS_NONE = 0, SS_BIP = 1, SS_BSQ = 2, };
    SAMP_STOR GetSampleStorage() {
        if (!m_pImgSet) return SS_NONE;
        const char* str = m_pImgSet->GetMetadataItem("INTERLEAVE", "IMAGE_STRUCTURE");
        if (!str) return SS_NONE;
        if (!strcmp(str, "PIXEL")) return SS_BIP;
        if (!strcmp(str, "BAND")) return SS_BSQ;
        return SS_NONE;
    }

    bool IsLoaded() { return m_pImgSet ? true : false; }
    const char* GetImagePath() { return m_strImagePath; }

    bool Open(const char* lpstrPath, GDALAccess openFlags)//,UINT flag = CImageBase::modeRead
    {
//      Close();
//      m_openFlags = (flag == CImageBase::modeRead) ? GA_ReadOnly : GA_Update;
        m_openFlags = openFlags;//GA_ReadOnly;
        m_pImgSet = (GDALDataset *)GDALOpen(lpstrPath, m_openFlags);
        if (!m_pImgSet) {
//          LogPrint(ERR_FILE_FLAG, GetLastErrorMsg()); 
            return false;
        }
        m_dataType = GetRasterBand(0)->GetRasterDataType();

//      if (m_dataType == GDT_UInt16 || m_dataType == GDT_Int16 ){
//          int nBands = GetBands();
//          m_tab16to8 = new int[nBands][65536];
//          if (m_tab16to8) {
//              for (int i = 0; i < nBands; i++){
//                  Get16to8Tab(m_tab16to8[i], i);
//              }
//          }
//      }

        if (!ReadTfw(lpstrPath)) {
            double adfGeoTransform[6];
            if (m_pImgSet->GetGeoTransform(adfGeoTransform) != CE_Failure)
            {
                UpdateThisGeoTransform(adfGeoTransform,false);
            }               
        }
        SetImagePath(lpstrPath);

        m_16to8_form = PERCENT;
        m_16to8_form_last = PERCENT;
        GetDefaultParam(m_16to8_form, &m_16to8_par1, &m_16to8_par2);

        return true;
    }
    bool Create(const char* lpstrPath, int nCols, int nRows, int nBands, GDALDataType dataType, SAMP_STOR samp_stor = SS_BIP)
    {
//      Close();
        GDALDataType defaultType;
        const char* pExt = strrchr(lpstrPath, '.'); pExt = GetDescription(pExt, &defaultType);
        if (!pExt) { LogPrint(ERR_FILE_FLAG, "unsupported file type %s", lpstrPath); return false; }
        if (dataType == GDT_Unknown) dataType = defaultType;
        GDALDriver *poDriver;
        poDriver = GetGDALDriverManager()->GetDriverByName(pExt);
        char **papszOptions = NULL;
        if (IsTiffDrive(pExt)){
            switch (samp_stor)
            {
            case image::SS_NONE:
                break;
            case image::SS_BIP:
                papszOptions = CSLSetNameValue(papszOptions, "INTERLEAVE", "PIXEL");
                break;
            case image::SS_BSQ:
                papszOptions = CSLSetNameValue(papszOptions, "INTERLEAVE", "BAND");
                break;
            default:
                break;
            }
            if (nBands>=3){
                papszOptions = CSLSetNameValue(papszOptions, "PHOTOMETRIC", "RGB");
            }
        }

        char **papszMetadata = poDriver->GetMetadata();
        if (CSLFetchBoolean(papszMetadata, GDAL_DCAP_CREATE, FALSE)){
            m_pImgSet = poDriver->Create(lpstrPath, nCols, nRows, nBands, dataType, papszOptions);
        }
        else{
//          char strPath[512];  strcpy(strPath, lpstrPath); strcpy(strrchr(strPath, '.'), ".tmp");
            m_pMemDriver = GetGDALDriverManager()->GetDriverByName("MEM");
            if (m_pMemDriver == NULL) { return false; }
            m_pImgSet = m_pMemDriver->Create("", nCols, nRows, nBands, dataType, papszOptions);//strPath
        }
        
        
        CSLDestroy(papszOptions);
        if (!m_pImgSet) return false;       

        m_dataType = dataType;
        m_openFlags = GA_Update;
        SetImagePath(lpstrPath);

        if (nBands >= 3){
            GetRasterBand(0)->SetColorInterpretation(GCI_RedBand);
            GetRasterBand(1)->SetColorInterpretation(GCI_GreenBand);
            GetRasterBand(2)->SetColorInterpretation(GCI_BlueBand);
//          if (nBands == 4){
//              GetRasterBand(3)->SetColorInterpretation(GCI_AlphaBand);
//          }
        }

        return true;
    }
    void Close()
    {
        Reset();
    }
    int GetRows() { return m_pImgSet->GetRasterYSize(); }
    int GetCols() { return m_pImgSet->GetRasterXSize(); }
    int GetBands(){ return m_pImgSet->GetRasterCount(); }
    int GetBits(int band = 0) {
        GDALDataType dataType = GetRasterBand(band)->GetRasterDataType();
        return GDALGetDataTypeSize(dataType);
    }
    int GetByte(int band = 0){
        return GetBits(band) / 8;
    }
    GDALDataType    GetDataType(){
        return m_dataType;
    }
    bool    GetStatistics(int nBandIdx, double *pdfMin, double *pdfMax, double *pdfMean, double *padfStdDev){
        if( GetRasterBand(nBandIdx)->GetStatistics(FALSE, TRUE, pdfMin, pdfMax, pdfMean, padfStdDev) == CE_None){
            return true;
        }
        return false;       
    }
    bool    GetNoDataValue(int nBandIdx, double* val){
        int bSuccess = 0;
        *val = GetRasterBand(nBandIdx)->GetNoDataValue(&bSuccess);
        return bSuccess ? true : false;
    }
    bool    SetNoDataValue(int nBandIdx, double val){
        if ( GetRasterBand(nBandIdx)->SetNoDataValue(val) == CE_None ){
            return true;
        }
        return false;
    }

    bool ImageIO(GDALRWFlag flag, void* data, int stCol, int stRow, int nCols, int nRows,
        int *panBandMap = NULL,
        int nPixelSpace = 0, int nLineSpace = 0, int nBandSpace = 0,
        float zoom = 1.0f, int* zoomCols = NULL, int* zoomRows = NULL)
    {
        int nCols_zoom = int(nCols*zoom);   int nRows_zoom = int(nRows*zoom);
        int nBytes = GetBits() / 8;
        if (!RasterIO(flag, stCol, stRow, nCols, nRows, data, nCols_zoom, nRows_zoom, m_dataType, GetBands(), panBandMap, nPixelSpace*nBytes, nLineSpace*nBytes, nBandSpace*nBytes))
            return false;
        if (zoomCols) *zoomCols = nCols_zoom;   if (zoomRows) *zoomRows = nRows_zoom;
        return true;
    }
    bool BandIO(GDALRWFlag flag, void* data, int band, int stCol, int stRow, int nCols, int nRows, float zoom = 1.0f, int* zoomCols = NULL, int* zoomRows = NULL)
    {
        int nCols_zoom = int(nCols*zoom);   int nRows_zoom = int(nRows*zoom);
        if ( !RasterBandIO(flag,band, stCol, stRow, nCols, nRows, data, nCols_zoom, nRows_zoom, m_dataType, 0, 0) )
            return false;
        if (zoomCols) *zoomCols = nCols_zoom;   if (zoomRows) *zoomRows = nRows_zoom;
        return true;
    }
    bool ImageIO_SampStore(GDALRWFlag flag, void* data, int stCol, int stRow, int nCols, int nRows, SAMP_STOR samp_store, bool bInvert3Bands)
    {
        int nBands = GetBands();    if (nBands < 2) return ImageIO(flag,data, stCol, stRow, nCols, nRows);
        
        if (samp_store==SS_NONE){
            return ImageIO(flag, data, stCol, stRow, nCols, nRows);
        }
        else {
            int *panBandMap = NULL, nPixelSpace, nLineSpace, nBandSpace;
            if (samp_store == SS_BSQ){
                if (!bInvert3Bands) return ImageIO(flag, data, stCol, stRow, nCols, nRows);
                nPixelSpace = 1;
                nLineSpace = nCols;
                nBandSpace = nCols*nRows;
            }
            else{
                nPixelSpace = nBands;
                nLineSpace = nCols*nPixelSpace;
                nBandSpace = 1;
            }
            panBandMap = new int[nBands];   for (int i = 0; i < nBands; i++) { panBandMap[i] = i + 1; }
            if (bInvert3Bands&&nBands >= 3) { panBandMap[0] = 3; panBandMap[1] = 2; panBandMap[2] = 1; }
            bool ret = ImageIO(flag, data, stCol, stRow, nCols, nRows, panBandMap, nPixelSpace, nLineSpace, nBandSpace);
            delete[] panBandMap;
            return ret;
        }
        return false;
    }
    bool BuildOverviews(const char *pszResampling,
        int nReqOverviews, int *panOverviewList,
        GDALProgressFunc pfnProgress,
        void *pProgressData)
    {
//      BuildOverviews(_T("NEAREST"), nLevel, pBandList, 0, NULL, GdalBuildPyramidProgress, NULL);
        return false;
    }
    
    bool        IsValid_16to8Tab() { return m_tab16to8 ? true : false; }
    BYTE        Tab_16to8(int val,int band = 0){ 
        return (BYTE)m_tab16to8[band][val];
    }

    bool Get16to8Tab(int nTab[65536],int band)
    {
        return Get16to8Tab(nTab, band, m_16to8_form, m_16to8_par1, m_16to8_par2);
//      int loSum, hiSum, sum = 0;
//      memset(nTab, 0, sizeof(int)* 65536);
//      if( !GetHistogram(nTab, band) ) return false;
//      int i, minI, maxI, midI, minSum, maxSum;
//      for (i = 0; i < 65535; i++) { sum += nTab[i]; }
// 
//      minSum = int(sum*0.001f);  if (minSum<10) minSum = 10;
//      maxSum = int(sum*0.001f);  if (maxSum<10) maxSum = 10;
//      for (sum = 0, minI = 0; minI<65535; minI++){ sum += nTab[minI]; if (sum>minSum) break; }
//      for (sum = 0, maxI = 65535; maxI>0; maxI--){ sum += nTab[maxI]; if (sum>maxSum) break; }
// 
// // for Sta. White rgn (Mybe Cloud)
// //                         float per = 1; midI = (minI + maxI) / 2;
// //                         for (loSum = 0, i = minI; i <= midI; i++) loSum += pTab16[i];
// //                         for (hiSum = 0, i = midI; i < maxI; i++) hiSum += pTab16[i];
// //                         per = (loSum*1.f / hiSum*1.f);
// //                         if (per < 0.05f || (maxI - minI) < 192 || maxI<160 || minI>768) ret = 10;
// 
//      memset(nTab, 0, sizeof(int)* 65536);
// //       printf("[band=%d] minI = %d\tmaxI = %d\n", band, minI, maxI);
//      if (minI == maxI){
//          nTab[minI] = maxI >> 8;
//      }
//      else{
//          for (i = minI; i < maxI; i++) nTab[i] = 1 + int((i - minI)*254.f / (maxI - minI));
//          for (; i < 65535; i++) nTab[i] = 255;
//      }
//      
// 
//      return true;
    }
    //对数,指数,百分比截取,标准截取,直方图均衡
    enum POPETATORFORM{LOGARITHM = 0,INDEX,PERCENT,BALANCE,STDPERCENT,WAVELET};
    void    ParseParam(int argc,const char **argv)
    {
        const char*  str16to8_form[] = { "log", "exp", "per", "bal", "stdp","wav" };
        if (m_tab16to8) delete[] m_tab16to8; m_tab16to8 = NULL;
        for (int i = 0; i < argc; i++){
            if (!stricmp(argv[i], "-pf")){
                if (++i >= argc) break;
                for (int j = 0; j < sizeof(str16to8_form) / sizeof(str16to8_form[0]); j++){
                    if (!stricmp(argv[i], str16to8_form[j])){                       
                        m_16to8_form =(POPETATORFORM) j;
                        GetDefaultParam(m_16to8_form, &m_16to8_par1, &m_16to8_par2);
                    }
                }           
            }
            else if (!stricmp(argv[i], "-par1")){
                if (++i >= argc) break;
                double a = atof(argv[i]);   if (a<0) continue;
                m_16to8_par1 = a;
            }
            else if (!stricmp(argv[i], "-par2")){
                if (++i >= argc) break;
                double a = atof(argv[i]);   if (a < 0) continue;
                m_16to8_par2 = a;
            }
        }
    }
    void        GetDefaultParam(POPETATORFORM form, double* par1, double* par2)
    {
        switch (form)
        {
        case LOGARITHM:
            *par1 = -1;     *par2 = -1;
            break;
        case INDEX:
            *par1 = 1;          *par2 = -1;
            break;
        case PERCENT:
            *par1 = 0.0025; *par2 = 0.0025;
            break;
        case BALANCE:
            *par1 = -1;         *par2 = -1;
            break;
        case STDPERCENT:
            *par1 = -1;         *par2 = -1;
            break;
        case WAVELET:
            *par1 = -1;         *par2 = -1;
            break;
        default:
            break;
        }
    }
    //灰度压缩方式,波段,参数1,参数2
    /************************************************************************/
    /* 参数1:log--无,ind--指数,per--头位百分比,bal--无,stdp--无               */
    /* 参数2:log--无,ind--无,per--尾位百分比,bal--无 ,stdp--无                  */
    /************************************************************************/
    bool Get16to8Tab(int nTab[65536], int band, POPETATORFORM form, double Par1, double Par2 = 0){  
        if (Par1 < 0 || Par2 < 0) Par1 = abs(Par1), Par2 = abs(Par2);
        const int bit16 = 65536, bit8 = 256;
        int i, j;
        double Imax, Imin, Imean, Istd;
        if (GetRasterBand(band)->ComputeStatistics(true, &Imin, &Imax, &Imean, &Istd, NULL, NULL) != CE_None) return false;
        for (i = 0; i < Imin; i++) nTab[i] = 0;
        for (i = Imax; i < bit16; i++) nTab[i] = 255;
        double nor, miners;
        switch (form)
        {
        case LOGARITHM:
            nor = 255.0 / log(Imax - Imin);
            nTab[static_cast<int>(Imin)] = 0;
            for (i = Imin + 1; i < Imax; i++)               nTab[i] = log(i - Imin) *nor;
            break;
        case INDEX:
            nor = 255.0 / pow(Imax - Imin, Par1);
            for (i = Imin; i < Imax; i++)
                nTab[i] = pow(i - Imin, Par1)*nor;
            break;
        case PERCENT:
            if (!GetHistogram(nTab, band))
                return false;
            {
                int loSum, hiSum, sum = 0;
                int i, minI, maxI, midI, minSum, maxSum;
                if (nTab[0]>1000) nTab[0] = 0;

                if (Par1 < 1 && Par2 < 1){
                    for (i = 0; i < 65535; i++) { sum += nTab[i]; }
                    minSum = int(sum*Par1);  if (minSum<10) minSum = 10;
                    maxSum = int(sum*Par2);  if (maxSum<10) maxSum = 10;
                    for (sum = 0, minI = 0; minI<65535; minI++){ sum += nTab[minI]; if (sum>minSum) break; }
                    for (sum = 0, maxI = 65535; maxI>0; maxI--){ sum += nTab[maxI]; if (sum>maxSum) break; }
                }
                else{
                    minI = int(Par1);
                    maxI = int(Par2);
//                  if (Imax>255 && Par2>Imax){
//                      maxI = Imax;
//                  }
                }
                

                memset(nTab, 0, sizeof(int)* 65536);
                if (minI == maxI){
                    nTab[minI] = maxI >> 8;
                }
                else{
                    for (i = minI; i < maxI; i++) nTab[i] = 1 + int((i - minI)*254.f / (maxI - minI));
                    for (; i < 65535; i++) nTab[i] = 255;
                }
            }
            
//          double *probability;//归一化概率分布
//          probability = new double[static_cast<int>(Imax - Imin)];
//          nor = 0;
//          for (i = Imin; i < Imax; i++) nor += nTab[i];
//          for (i = Imin; i < Imax; i++) probability[i - int(Imin)] = 1.0*nTab[i] / nor;
//          miners = 0;
//          //首截取
//          for (i = Imin;; i++){
//              if (miners > Par1)
//                  break;
//              miners += probability[i - int(Imin)];
//              nTab[i] = 0;
//          }
//          miners = 0;
//          //尾截取
//          for (j = Imax;; j--){
//              if (miners > Par2)
//                  break;
//              miners += probability[j - int(Imin)];
//              nTab[j] = 255;
//          }
//          //线性拉伸
//          for (int k = i; k < j + 1; k++)
//              nTab[k] = 255 * (k - i) / (j - i);
//          for (i = Imax; i < bit16; i++) nTab[i] = 255;
//          i = 0;
            break;
        case STDPERCENT:
            if (Imin < Imean - 3 * Istd)
                miners = Imean - 3 * Istd;
            else miners = Imin;
            nor = 3 * Istd + Imean - miners;
            nor = 255.0 / nor;
            for (i = Imin; i < Imean - 3 * Istd; i++)               nTab[i] = 0;
            for (; i < Imean + 3 * Istd; i++)               nTab[i] = (i - miners)*nor;
            for (; i < Imax; i++)               nTab[i] = 255;
            break;
        case BALANCE:
            if (!GetHistogram(nTab, band))
                return false;
            nor = 0;
            for (i = Imin; i < Imax; i++) nor += nTab[i];
            miners = 0;
            for (i = Imin; i < Imax; i++){
                miners += nTab[i] / nor;
                nTab[i] = miners * 255;
            }
            for (i = Imax; i < bit16; i++) nTab[i] = 255;
            break;
        case WAVELET:
            int waveletSacle[] = { 1024, 512, 256 };
            if (Imin < Imean - 3 * Istd)
                miners = Imean - 3 * Istd;
            else miners = Imin;
            nor = 3 * Istd + Imean - miners;
            if (!GetHistogram(nTab, band))
                return false;
            double nTabTemp[2048] = { 0 }, detail[2048];
            int dstTab[65536] = { 0 };//目标色表
            unsigned int nTabSum = 0;
            //ntab赋值
            for (i = Imin; i < Imean - 3 * Istd; i++)               nTab[i] = 0;
            for (; i < Imean + 3 * Istd; i++)               nTabSum += nTab[i];
            for (; i < Imax; i++)               nTab[i] = 255;
            //待压缩数组赋值(拉伸至1024)
            for (i = miners; i < Imean + 3 * Istd; i++)
                nTabTemp[static_cast<int>((i - miners) / nor * 1024)] = 1.0*nTab[i] / nTabSum;
            //db3小波变换
            DWT(6, nTabTemp, detail, 3, waveletSacle);
            nor = 0;
            //目标区域归一化
            for (i = 1536; i < 1792; i++)       nor += nTabTemp[i];
            for (i = 1536; i < 1792; i++)       nTabTemp[i] /= nor;
            //直方图匹配
            double accO = 0, accT = 0;
            for (i = 1536, j = miners; accT < 1 - 1e-9 || accO < 1 - 1e-9;){
                dstTab[j] = i - 1536;
                if (accO < accT)
                    accO += 1.0*nTab[j] / nTabSum, j++;//目标直方图
                else accT += nTabTemp[i], i++;//原始直方图
            }
            while (j != 65536){
                dstTab[j] = 255;
                j++;
            }
            memcpy(nTab, dstTab, 65536);
            break;
        }
        return true;
    }
    bool Read8(BYTE* data, int stCol, int stRow, int nCols, int nRows, 
        int *panBandMap = NULL,
        int nPixelSpace = 0, int nLineSpace = 0, int nBandSpace = 0, 
        float zoom = 1.0f, int* zoomCols = NULL, int* zoomRows = NULL)
    {
        switch (m_dataType)
        {
        case  GDT_Byte:
            return ImageIO(GF_Read,data, stCol, stRow, nCols, nRows, panBandMap, nPixelSpace, nLineSpace, nBandSpace, zoom, zoomCols, zoomRows);
        case GDT_UInt16:
        case GDT_Int16:
        {
                          if (!m_tab16to8) {
                              int nBands = GetBands();
                              m_tab16to8 = new int[nBands][65536];
                              for (int i = 0; i < nBands; i++){
                                  Get16to8Tab(m_tab16to8[i], i);
                              }
                          }
                          int nCols_new = int(nCols*zoom);
                          int nRows_new = int(nRows*zoom);
                          int rowlen = nCols_new*GetBands();
                          int bandSpace = nCols_new*nRows_new;
                          int nBytes = 2;
                          WORD* pRow = new WORD[rowlen + 8];
                          WORD* pT; int r, c, b;

                          for (r = 0; r < nRows_new; r++)
                          {   
                              if (!RasterIO(GF_Read, stCol, stRow + int(r / zoom + 0.5), nCols, 1, pRow, nCols_new, 1, m_dataType, GetBands(), panBandMap, nPixelSpace*nBytes, nLineSpace*nBytes, nBandSpace*nBytes))
                              {
                                  if (!panBandMap){
                                      for (b = 0; b < GetBands(); b++) memset(data + bandSpace*b, 0, sizeof(BYTE)*nCols_new);
                                      data += nCols_new;
                                  }
                                  else{
                                      memset(data, 0, rowlen);
                                      data += rowlen;
                                  }
                                  continue;
                              }
                              if (!panBandMap){
                                  for (b = 0; b < GetBands(); b++){
                                      BYTE* pData = data + bandSpace*b;
                                      pT = pRow + b*nCols_new;
                                      for (c = 0; c < nCols_new; c++, pT++, pData++){ *pData = Tab_16to8(*pT,b); }
                                  }
                                  data += nCols_new;
                              }
                              else{//error:Tab_16to8(band)
                                  for (c = 0, pT = pRow; c < rowlen; c++, pT++, data++){ *data = Tab_16to8(*pT); }
                              }                           
                          }
                          delete[] pRow;
//                        delete[] pTab16;
                          if (zoomCols) *zoomCols = nCols_new;  if (zoomRows) *zoomRows = nRows_new;
                          return true;
        }
        default:
            return false;
        }
    }
    bool ReadBand8(BYTE* data, int band, int stCol, int stRow, int nCols, int nRows, float zoom = 1.0f, int* zoomCols = NULL, int* zoomRows = NULL)
    {
        switch (m_dataType)
        {
        case  GDT_Byte:
            return BandIO(GF_Read,data,band, stCol, stRow, nCols, nRows, zoom,zoomCols,zoomRows);
        case GDT_UInt16:
        case GDT_Int16:
        {
                          if (!m_tab16to8) {
                              int nBands = GetBands();
                              m_tab16to8 = new int[nBands][65536];
                              for (int i = 0; i < nBands; i++){
                                  Get16to8Tab(m_tab16to8[i], i);
                              }
                          }
                          int nCols_new = int(nCols*zoom);
                          int nRows_new = int(nRows*zoom);
                          WORD* pRow = new WORD[nCols_new + 8];
                          WORD* pT; int r, c;
                          
                          for ( r = 0; r < nRows_new; r++)
                          {
                              if (!RasterBandIO(GF_Read, band, stCol, stRow + int(r / zoom + 0.5), nCols, 1, pRow, nCols_new, 1, m_dataType, 0, 0))
                              {
                                  memset(data, 0, sizeof(BYTE)*nCols_new);  data += nCols_new;
                                  continue;
                              }
                               
                              for (c = 0, pT = pRow; c < nCols_new; c++, pT++, data++){ 
                                  *data = Tab_16to8(*pT,band); 
                              }
                          }
                          delete[] pRow;
//                        delete[] pTab16;
                          if (zoomCols) *zoomCols = nCols_new;  if (zoomRows) *zoomRows = nRows_new;
                          return true;
        }
        default:
            return false;
        }
    }
    bool ReadGray8(BYTE* data, int stCol, int stRow, int nCols, int nRows, float zoom = 1.0f, int* zoomCols = NULL, int* zoomRows = NULL,bool bBGR = true)
    {
        if (GetBands() == 1){
//          return Read8(data, stCol, stRow, nCols, nRows, 0, 0, 0, 0, zoom, zoomCols, zoomRows);
            return ReadBand8(data, 0, stCol, stRow, nCols, nRows, zoom, zoomCols, zoomRows);
        }
        else if (GetBands() >= 3){
            int nCols_new = int(nCols*zoom);
            int nRows_new = int(nRows*zoom);
            GDALRasterBand *pBandB, *pBandG, *pBandR;
            int bandB, bandG, bandR;
            if (bBGR) {             
                bandB = 0;  bandG = 1;  bandR = 2;
            }
            else
            {
                bandB = 2;  bandG = 1;  bandR = 0;
            }
            pBandB = GetRasterBand(bandB);
            pBandG = GetRasterBand(bandG);
            pBandR = GetRasterBand(bandR);
            switch (m_dataType)
            {
            case  GDT_Byte:
            {                             
                              BYTE* pRowB = new BYTE[nCols_new + 8];
                              BYTE* pRowG = new BYTE[nCols_new + 8];
                              BYTE* pRowR = new BYTE[nCols_new + 8];
                              int r, c; BYTE *pB, *pR, *pG;
                              for (r = 0; r < nRows_new; r++)
                              {
                                  if (pBandB->RasterIO(GF_Read, stCol, stRow + int(r / zoom + 0.5), nCols, 1, pRowB, nCols_new, 1, m_dataType, 0, 0)==CE_Failure)
                                  {
                                      memset(pRowB, 0, sizeof(BYTE)*nCols_new); 
                                      continue;
                                  }
                                  if (pBandG->RasterIO(GF_Read, stCol, stRow + int(r / zoom + 0.5), nCols, 1, pRowG, nCols_new, 1, m_dataType, 0, 0) == CE_Failure)
                                  {
                                      memset(pRowG, 0, sizeof(BYTE)*nCols_new);
                                      continue;
                                  }
                                  if (pBandR->RasterIO(GF_Read, stCol, stRow + int(r / zoom + 0.5), nCols, 1, pRowR, nCols_new, 1, m_dataType, 0, 0) == CE_Failure)
                                  {
                                      memset(pRowR, 0, sizeof(BYTE)*nCols_new);
                                      continue;
                                  }
                                  for (c = 0,pB = pRowB,pG = pRowG,pR = pRowR; c < nCols_new; c++,pB++,pG++,pR++, data++){ 
                                      *data = BGR2GRAY(*pB, *pG, *pR);
                                  }
                              }
                              delete[] pRowG;
                              delete[] pRowB;
                              delete[] pRowR;
                              if (zoomCols) *zoomCols = nCols_new;  if (zoomRows) *zoomRows = nRows_new;
                              return true;
            }
            case GDT_UInt16:
            case GDT_Int16:
            {
                              if (!m_tab16to8) {
                                  int nBands = GetBands();
                                  m_tab16to8 = new int[nBands][65536];
                                  for (int i = 0; i < nBands; i++){
                                      Get16to8Tab(m_tab16to8[i], i);
                                  }
                              }
                              WORD* pRowB = new WORD[nCols_new + 8];
                              WORD* pRowG = new WORD[nCols_new + 8];
                              WORD* pRowR = new WORD[nCols_new + 8];
                              int r, c; WORD *pB, *pR, *pG;
                              for (r = 0; r < nRows_new; r++)
                              {
                                  if (pBandB->RasterIO(GF_Read, stCol, stRow + int(r / zoom + 0.5), nCols, 1, pRowB, nCols_new, 1, m_dataType, 0, 0) == CE_Failure)
                                  {
                                      memset(pRowB, 0, sizeof(WORD)*nCols_new);
                                      continue;
                                  }
                                  if (pBandG->RasterIO(GF_Read, stCol, stRow + int(r / zoom + 0.5), nCols, 1, pRowG, nCols_new, 1, m_dataType, 0, 0) == CE_Failure)
                                  {
                                      memset(pRowG, 0, sizeof(WORD)*nCols_new);
                                      continue;
                                  }
                                  if (pBandR->RasterIO(GF_Read, stCol, stRow + int(r / zoom + 0.5), nCols, 1, pRowR, nCols_new, 1, m_dataType, 0, 0) == CE_Failure)
                                  {
                                      memset(pRowR, 0, sizeof(WORD)*nCols_new);
                                      continue;
                                  }
                                  for (c = 0, pB = pRowB, pG = pRowG, pR = pRowR; c < nCols_new; c++, pB++, pG++, pR++, data++){
                                      *data = BGR2GRAY(Tab_16to8(*pB, bandB), Tab_16to8(*pG, bandG), Tab_16to8(*pR,bandR));
                                  }
                              }
                              delete[] pRowG;
                              delete[] pRowB;
                              delete[] pRowR;
                              if (zoomCols) *zoomCols = nCols_new;  if (zoomRows) *zoomRows = nRows_new;
                              return true;
            }
            default:
                return false;
            }
        }
        return false;
    }

    bool GetGeoTransform(double * padfTransform) { memcpy(padfTransform, m_adfGeoTransform, sizeof(m_adfGeoTransform));  return true; }
    const char *GetProjectionRef(void) { return m_pImgSet->GetProjectionRef(); }

    bool SetGeoInformation(double * padfTransform,const char* lpstrProjectionRef,bool bPixelCenter = true) {
        bool ret = true;
        if (padfTransform)  {
            UpdateThisGeoTransform(padfTransform,bPixelCenter);
            if (bPixelCenter) {
                GEOTRANSFROM_CENTER2CORNER(padfTransform);
            }
        }
        if (lpstrProjectionRef){
            WriteProj(m_strImagePath, lpstrProjectionRef);
        }
        if (IsLoaded()) {
            if (m_openFlags == GA_Update)   {
                if (padfTransform)  {
                    if (m_pImgSet->SetGeoTransform(padfTransform) == CE_Failure) ret = false;
                    WriteTfw(m_strImagePath);
                }
                if (lpstrProjectionRef) 
                    if (m_pImgSet->SetProjection(lpstrProjectionRef) == CE_Failure) ret = false;
            }
            else
            {
                if (m_pImgSet) delete m_pImgSet;    m_pImgSet = NULL;
                m_pImgSet = (GDALDataset *)GDALOpen(m_strImagePath, GA_Update);
                if (padfTransform)  {
                    if (m_pImgSet->SetGeoTransform(padfTransform) == CE_Failure) ret = false;
                    WriteTfw(m_strImagePath);
                }
                if (lpstrProjectionRef) 
                    if (m_pImgSet->SetProjection(lpstrProjectionRef) == CE_Failure) ret = false;
                delete m_pImgSet;
                m_pImgSet = (GDALDataset *)GDALOpen(m_strImagePath, m_openFlags);
            }
        }
        return ret;
    }
//  bool UpdateGeoTransform(double * padfTransform)
//  {
//      memcpy(m_adfGeoTransform, padfTransform, sizeof(m_adfGeoTransform));
//  }
    void GetGeoXY(double* x, double* y)
    {
        GDALApplyGeoTransform(m_adfGeoTransform,
            *x, *y, x, y);
    }
    void GetImgXY(double* x, double* y)
    {
        GDALApplyGeoTransform(m_adfReverseGeoTransform,
            *x, *y, x, y);
    }

    const char* GetLastErrorMsg() { return CPLGetLastErrorMsg(); }
protected:
    const char* GetDescription(const char* lpstrExt,GDALDataType* defaultType)
    {
        static const char* strExt[] = { "bmp", "jpg", "tif", "tiff", "img", "bt", "ecw", "fits", "gif", "hdf", "hdr","pix","png" };
        static const char* strDescrip[] = { "BMP", "JPEG", "GTiff","GTiff", "HFA", "BT", "ECW", "FITS", "GIF", "HDF4", "EHdr","PCIDSK","PNG" };
        static GDALDataType dataType[] = { GDT_Byte, GDT_Byte, GDT_Byte, GDT_Byte, GDT_Byte, GDT_Byte, GDT_Byte, GDT_Byte, GDT_Byte, GDT_Byte, GDT_Byte };

        char suffix[10];    if (*lpstrExt == '.') strcpy(suffix, lpstrExt + 1); else strcpy(suffix, lpstrExt);
        strlwr(suffix);

        for (unsigned long i = 0; i < sizeof(strExt) / sizeof(const char*);i++){
            if (!strcmp(strExt[i], suffix)) {
                if (defaultType) *defaultType = dataType[i];
                return strDescrip[i]; 
            }
        }
        return NULL;
    }
    bool       IsTiffDrive(const char* descrip){
        if (!strcmp(descrip, "GTiff")) return true;
        return false;
    }
    void Reset()
    {
//      m_openFlags = -1;
        
        if (m_pMemDriver){
            if (strlen(m_strImagePath)>3 && m_pImgSet){
                GDALDriver *pDstDriver = NULL;
                GDALDataType defaultType;
                const char* pExt = strrchr(m_strImagePath, '.');    pExt = GetDescription(pExt, &defaultType);
                if (pExt){
                    pDstDriver = GetGDALDriverManager()->GetDriverByName(pExt);
                    if (pDstDriver){
                        GDALDataset* pDst = pDstDriver->CreateCopy(m_strImagePath, m_pImgSet, FALSE, NULL, NULL, NULL);
                        if (pDst) delete pDst;
                    }
                }       
            }
            m_pMemDriver = NULL;
        }
        memset(m_strImagePath, 0, 512);
        if (m_pImgSet) delete m_pImgSet;    m_pImgSet = NULL;
        if (m_tab16to8) delete[] m_tab16to8;    m_tab16to8 = NULL;
        m_adfGeoTransform[0] = 0.0;
        m_adfGeoTransform[1] = 1.0;
        m_adfGeoTransform[2] = 0.0;

        m_adfGeoTransform[3] = 0.0;
        m_adfGeoTransform[4] = 0.0;
        m_adfGeoTransform[5] = 1.0;

        memcpy(m_adfReverseGeoTransform, m_adfGeoTransform, sizeof(m_adfReverseGeoTransform));
    }
    bool WriteProj(const char* lpstrPathName, const char* lpstrWkt, const char* lpstrExt = "prj")
    {
        char strPath[512];  strcpy(strPath, lpstrPathName); strcpy(strrchr(strPath, '.') + 1, lpstrExt);
        FILE* fp = fopen(strPath, "w"); if (!fp) return false;
        fprintf(fp, "%s", lpstrWkt);
        fclose(fp);
        return true;
    }
    bool ReadTfw(const char* lpstrPathName,const char* lpstrExt = "tfw")
    {
        char strPath[512];  strcpy(strPath, lpstrPathName); strcpy(strrchr(strPath, '.') + 1, lpstrExt);
        FILE* fp = fopen(strPath, "rt");    if (!fp) return false;
        double adfGeoTransform[6];
        fscanf(fp, "%lf%lf%lf%lf%lf%lf", 
            adfGeoTransform + 1, adfGeoTransform + 2,
            adfGeoTransform + 4, adfGeoTransform + 5,
            adfGeoTransform + 0, adfGeoTransform + 3);
        fclose(fp);
        UpdateThisGeoTransform(adfGeoTransform);
        return true;
    }
    bool WriteTfw(const char* lpstrPathName, const char* lpstrExt = "tfw")
    {
        char strPath[512];  strcpy(strPath, lpstrPathName); strcpy(strrchr(strPath, '.') + 1, lpstrExt);
        FILE* fp = fopen(strPath, "wt");    if (!fp) return false;
        fprintf(fp, "\t\t%lf\n\t\t%lf\n\t\t%lf\n\t\t%lf\n\t\t%lf\n\t\t%lf\n",
            m_adfGeoTransform[1], m_adfGeoTransform[2],
            m_adfGeoTransform[4], m_adfGeoTransform[5],
            m_adfGeoTransform[0], m_adfGeoTransform[3]);
        fclose(fp);
        return true;
    }
    void SetImagePath(const char* lpstrPathName)
    {
        strcpy(m_strImagePath, lpstrPathName);
    }
    void UpdateThisGeoTransform(double* padfTransform,bool bPixelCenter = true )
    {
        memcpy(m_adfGeoTransform, padfTransform, sizeof(m_adfGeoTransform));
        if (!bPixelCenter) GEOTRANSFROM_CORNER2CENTER(m_adfGeoTransform);
        GDALInvGeoTransform(m_adfGeoTransform, m_adfReverseGeoTransform);
    }
public:
    bool RasterIO(GDALRWFlag flag,int nXOff, int nYOff, int nXSize, int nYSize,
        void * pData, int nBufXSize, int nBufYSize,
        GDALDataType eBufType,
        int nBandCount, int *panBandMap,
        int nPixelSpace, int nLineSpace, int nBandSpace)
    {
        if (m_pImgSet->RasterIO(flag, nXOff, nYOff, nXSize, nYSize, pData, nBufXSize, nBufYSize, eBufType, nBandCount, panBandMap, nPixelSpace, nLineSpace, nBandSpace) == CE_Failure)
            return false;
        return true;
    }
    inline GDALRasterBand *GetRasterBand(int band){
        return m_pImgSet->GetRasterBand(band + 1);
    }
    bool RasterBandIO(GDALRWFlag flag, int band,
        int nXOff, int nYOff, int nXSize, int nYSize,
        void * pData, int nBufXSize, int nBufYSize,
        GDALDataType eBufType,
        int nPixelSpace,
        int nLineSpace)
    {
        if (GetRasterBand(band)->RasterIO(flag, nXOff, nYOff, nXSize, nYSize, pData, nBufXSize, nBufYSize, eBufType, nPixelSpace, nLineSpace) == CE_Failure)
            return false;
        return true;
    }
    bool GetHistogram(int pTab[65536], int band)
    {
//      LogPrint(0, "Calc histogram...\n");
#if GDAL_VERSION_NUM >= 2000000 
        GUIntBig pTabTemp[65536];
#else
        int* pTabTemp = pTab;
#endif
        if (GetRasterBand(band)->GetHistogram(-0.5, 65536.5, 65536, pTabTemp, 0, 0, NULL, NULL) == CE_Failure)//GDALTermProgress
            return false;
#if GDAL_VERSION_NUM >= 2000000 
        for (int i = 0; i < 65536; i++)
            pTab[i] = static_cast<int>(pTabTemp[i]);
#endif
        return true;
    }
protected:
    enum IMAGECOOR { left_top, left_bottom } m_imgCoor;
    GDALDataset*    m_pImgSet;
    GDALDriver*     m_pMemDriver;
    GDALAccess      m_openFlags;

    GDALDataType    m_dataType;

    double      m_adfGeoTransform[6];
    double      m_adfReverseGeoTransform[6];

    char            m_strImagePath[512];
    int             (*m_tab16to8)[65536];

    POPETATORFORM   m_16to8_form;
    POPETATORFORM   m_16to8_form_last;
    double                      m_16to8_par1;
    double                      m_16to8_par2;
};
// This is the constructor of a class that has been exported.
// see imagebase.h for the class definition
#define MAX_BUFSIZE     (64*1024*1024)
CImageBase::CImageBase()
{
    m_pImgSet = new image;
    m_pImgBuf = NULL;
//  m_pImgBuf = NULL;
    Reset();
    return;
}

CImageBase::~CImageBase()
{
    if (m_pImgSet) delete m_pImgSet;    m_pImgSet = NULL;
    if (m_pImgBuf) delete[] m_pImgBuf;    m_pImgBuf = NULL;
}

void CImageBase::Reset()
{
    if (m_pImgBuf) delete[] m_pImgBuf; m_pImgBuf = NULL;
    m_stBufCol = -99999;
    m_stBufRow = +99999;
    m_nBufCols = m_nBufRows = 0;
}
bool CImageBase::IsLoaded()
{
    return m_pImgSet ? m_pImgSet->IsLoaded() : false;
}
bool CImageBase::Open(const char* lpstrPath, int mode /* = modeRead */ )
{
    Close();
    if (m_pImgSet->Open(lpstrPath,mode==modeRead?GA_ReadOnly:GA_Update)){
        calc_buf_size(&m_nBufCols,&m_nBufRows,m_pImgSet->GetByte());
        return true;
    }
    return false;
}
bool CImageBase::Create(const char* lpstrPath, int nCols, int nRows, int nBands, int nBits /* = 8 */, SAMP_STOR samp_store /* = SS_BIP */ )
{
    Close();
    GDALDataType dataType;
    
    image::SAMP_STOR stor = image::SS_BIP;
    switch (samp_store)
    {
    case CImageBase::SS_NONE:
        break;
    case CImageBase::SS_BIP:
        stor = image::SS_BIP;
        break;
    case CImageBase::SS_BSQ:
        stor = image::SS_BSQ;
        break;
    default:
        break;
    }
    switch (nBits)
    {
    case 8:
        dataType = GDT_Byte; break;
    case 16:
        dataType = GDT_UInt16; break;
    case 32:
        dataType = GDT_Float32; break;
    case 64:
        dataType = GDT_Float64; break;
    default:
        dataType = GDT_Unknown;
        break;
    }
    if (m_pImgSet->Create(lpstrPath, nCols, nRows, nBands, dataType, stor)){
//      calc_buf_size(&m_nBufCols, &m_nBufRows);
        return true;
    }
    return false;
}
void CImageBase::Close()
{
    Reset();
    m_pImgSet->Close();
}
void CImageBase::ParseParam(int argc,const char **argv)
{
    m_pImgSet->ParseParam(argc, argv);
    return;
}
bool CImageBase::Read(void* data, int stCol, int stRow, int nCols, int nRows, float zoom /* = 1.0f */, int* zoomCols /* = NULL */, int* zoomRows /* = NULL */)
{
    return m_pImgSet->ImageIO(GF_Read,data, stCol, stRow, nCols, nRows, 0, 0, 0, 0, zoom, zoomCols, zoomRows);
}
bool CImageBase::Write(void* data, int stCol, int stRow, int nCols, int nRows, float zoom /* = 1.0f */, int* zoomCols /* = NULL */, int* zoomRows /* = NULL */)
{
    return m_pImgSet->ImageIO(GF_Write, data, stCol, stRow, nCols, nRows, 0, 0, 0, 0, zoom, zoomCols, zoomRows);
}
bool CImageBase::ReadBand(void* data, int band, int stCol, int stRow, int nCols, int nRows, float zoom /* = 1.0f */, int* zoomCols /* = NULL */, int* zoomRows /* = NULL */)
{
    return m_pImgSet->BandIO(GF_Read,data, band, stCol, stRow, nCols, nRows, zoom, zoomCols, zoomRows);
}
bool CImageBase::WriteBand(void* data, int band, int stCol, int stRow, int nCols, int nRows, float zoom /* = 1.0f */, int* zoomCols /* = NULL */, int* zoomRows /* = NULL */)
{
    return m_pImgSet->BandIO(GF_Write,data, band, stCol, stRow, nCols, nRows, zoom, zoomCols, zoomRows);
}
bool CImageBase::Read8(BYTE* data, int stCol, int stRow, int nCols, int nRows, float zoom /* = 1.0f */, int* zoomCols /* = NULL */, int* zoomRows /* = NULL */)
{
    return m_pImgSet->Read8(data, stCol, stRow, nCols, nRows, 0, 0, 0, 0, zoom, zoomCols, zoomRows);
}
bool CImageBase::ReadBand8(BYTE* data, int band, int stCol, int stRow, int nCols, int nRows, float zoom /* = 1.0f */, int* zoomCols /* = NULL */, int* zoomRows /* = NULL */)
{
    return m_pImgSet->ReadBand8(data,band, stCol, stRow, nCols, nRows, zoom,zoomCols,zoomRows);
}

bool CImageBase::ReadGray8(BYTE* data, int stCol, int stRow, int nCols, int nRows, float zoom /* = 1.0f */, int* zoomCols /* = NULL */, int* zoomRows /* = NULL */, bool bBGR /* = true */)
{
    return m_pImgSet->ReadGray8(data, stCol, stRow, nCols, nRows, zoom, zoomCols, zoomRows, bBGR);
}
bool CImageBase::Read(void* data, int stCol, int stRow, int nCols, int nRows, SAMP_STOR samp_store, bool bInvert3Bands)
{
    return m_pImgSet->ImageIO_SampStore(GF_Read, data, stCol, stRow, nCols, nRows, (image::SAMP_STOR)samp_store, bInvert3Bands);
}
bool CImageBase::Write(void* data, int stCol, int stRow, int nCols, int nRows, SAMP_STOR samp_store, bool bInvert3Bands)
{
    return m_pImgSet->ImageIO_SampStore(GF_Write, data, stCol, stRow, nCols, nRows, (image::SAMP_STOR)samp_store, bInvert3Bands);

}

int CImageBase::GetRows() const{ return m_pImgSet->GetRows(); }
int CImageBase::GetCols() const{ return m_pImgSet->GetCols(); }
int CImageBase::GetBands() const{ return m_pImgSet->GetBands(); }
int CImageBase::GetBits() const{ return m_pImgSet->GetBits(); }
int CImageBase::GetByte() const{ return m_pImgSet->GetByte(); }

void CImageBase::GetGeoTransform(double * padfTransform, bool bPixelCenter /* = true */) const
{
    if (m_pImgSet->GetGeoTransform(padfTransform)){
        if (!bPixelCenter) GEOTRANSFROM_CENTER2CORNER(padfTransform);
    }
}
const char* CImageBase::GetProjectionRef() const{ return m_pImgSet->GetProjectionRef(); }

bool CImageBase::SetGeoTransform(double * padfTransform,bool bPixelCenter /* = true */)
{
    return m_pImgSet->SetGeoInformation(padfTransform, NULL, bPixelCenter);
}
bool CImageBase::SetProjectionRef(const char* lpstrProjectionRef) { return m_pImgSet->SetGeoInformation(NULL, lpstrProjectionRef); }

#define EQU_ZERO(a) (fabs(a) <1e-6)
bool CImageBase::CopyGeoInformation(const CImageBase& img)
{
    double padfTransform[6] = {0,1,0,0,0,1};    img.GetGeoTransform(padfTransform);
    if (EQU_ZERO(padfTransform[0]) && EQU_ZERO(padfTransform[1] - 1) && EQU_ZERO(padfTransform[2]) && EQU_ZERO(padfTransform[3]) && EQU_ZERO(padfTransform[4]) && EQU_ZERO(padfTransform[5] - 1)) return false;
    return m_pImgSet->SetGeoInformation(padfTransform,img.GetProjectionRef());
}

const char* CImageBase::GetImagePath() const { return m_pImgSet->GetImagePath(); }

void   CImageBase::calc_buf_size(int* nCols, int* nRows,int datasize)
{
    int rows, cols;
    int bufsz_band = MAX_BUFSIZE / GetBands();

    rows = bufsz_band / (datasize*GetCols());   if (rows < 10) rows = 10;
    cols = bufsz_band / (datasize*rows);

    if (cols > GetCols()) cols = GetCols(); if (rows > GetRows()) rows = GetRows();

    *nCols = cols;  *nRows = rows;
}
bool    CImageBase::malloc_data_buf(int nCols, int nRows, int nBands, int datasize){
    if (m_pImgBuf) delete[] m_pImgBuf; m_pImgBuf = NULL;
    m_pImgBuf = new BYTE[nCols*nRows*nBands*datasize];
    if (!m_pImgBuf) return false;
    m_nBufRows = nRows; m_nBufCols = nCols;
    return true;
}

bool    CImageBase::adjust_data_buf(double fcol, double frow)
{   
    if (fcol >= m_stBufCol && frow >= m_stBufRow && fcol < m_stBufCol + m_nBufCols-1 && frow < m_stBufRow + m_nBufRows-1 ){
        return false;
    }
    int col, row;
    col = int(fcol);    row = int(frow);

    int half_cols = m_nBufCols / 2; int half_rows = m_nBufRows / 2;
    if (col - half_cols < 0) m_stBufCol = 0; else m_stBufCol = col - half_cols;
    if (row - half_rows < 0) m_stBufRow = 0; else m_stBufRow = row - half_rows;

    if (m_stBufCol + m_nBufCols>GetCols()) m_stBufCol = GetCols() - m_nBufCols;
    if (m_stBufRow + m_nBufRows>GetRows()) m_stBufRow = GetRows() - m_nBufRows;

    m_pImgSet->ImageIO(GF_Read, m_pImgBuf, m_stBufCol, m_stBufRow, m_nBufCols, m_nBufRows);
//  m_pImgSet->RasterIO(GF_Read, m_stBufCol, m_stBufRow, m_nBufCols, m_nBufRows, m_pImgBuf, m_nBufCols, m_nBufRows, eBufType, GetBands(), 0, 0, 0, 0);

    return true;
}

BYTE* CImageBase::GetBandBuf(int col, int row, int band){
    int datasize = m_pImgSet->GetByte();
    return m_pImgBuf + (band*m_nBufRows*m_nBufCols + row*m_nBufCols + col)*datasize;
}

void    CImageBase::SetBandBufVal(int col, int row, int band, double val){
#define SET_VALUE(type,buf,val) \
    { \
    *((type*)buf) = (type)val; \
    }\

    BYTE* pBuf = GetBandBuf(col, row, band);
    switch (m_pImgSet->GetDataType())
    {
    case GDT_UInt16:
        SET_VALUE(unsigned short,pBuf,val);
        break;
    case GDT_Int16:
        SET_VALUE(short, pBuf, val);
        break;
    case GDT_UInt32:
        SET_VALUE(unsigned int, pBuf, val);
        break;
    case GDT_Int32:
        SET_VALUE(int, pBuf, val);
        break;
    case GDT_Float32:
        SET_VALUE(float, pBuf, val);
        break;
    case GDT_Float64:
        SET_VALUE(double, pBuf, val);
        break;
    default:
        SET_VALUE(BYTE, pBuf, val);
        break;
    }
}

double  CImageBase::GetBandBufVal(double col, double row, int band)
{
#define BILINEAR_INTERPOLATION(type)        \
    { \
    type* pDataY0 = (type*)pBuf;    \
    type* pDataY1 = pDataY0 + m_nBufCols;   \
    if (dfDeltaX < 1e-5) {\
    if (dfDeltaY < 1e-5) return (double)*pDataY0; \
    return pDataY0[0] * dfDeltaY1 + pDataY1[0] * dfDeltaY;\
    }\
    if (dfDeltaY < 1e-5){ return pDataY0[0] * dfDeltaX1 + pDataY0[1] * dfDeltaX; }\
    if (pDataY0[0] == NOVALUE || pDataY0[1] == NOVALUE || pDataY1[0] == NOVALUE || pDataY1[1] == NOVALUE) \
    { if (pDataY0[0] != NOVALUE) return pDataY0[0]; if (pDataY0[1] != NOVALUE) return pDataY0[1]; \
    if (pDataY1[0] != NOVALUE) return pDataY1[0]; if (pDataY1[1] != NOVALUE) return pDataY1[1];  return NOVALUE; }  \
    dfXZ1 = pDataY0[0] * dfDeltaX1 + pDataY0[1] * dfDeltaX; \
    dfXZ2 = pDataY1[0] * dfDeltaX1 + pDataY1[1] * dfDeltaX; \
    return (dfXZ1 * dfDeltaY1 + dfXZ2 * dfDeltaY);  \
    }   \

    int nCol = int(col);    int nRow = int(row);

    double dfDeltaX = col - nCol;
    double dfDeltaX1 = 1.0 - dfDeltaX;
    double dfDeltaY = row - nRow;
    double dfDeltaY1 = 1.0 - dfDeltaY;

    BYTE* pBuf = GetBandBuf(nCol, nRow, band);

    double dfXZ1, dfXZ2;
    switch (m_pImgSet->GetDataType())
    {
//  case 1:
//      dfXZ1 = pBuf[0] * dfDeltaX1 + pBuf[1] * dfDeltaX;
//      pBuf += m_nBufCols;
//      dfXZ2 = pBuf[0] * dfDeltaX1 + pBuf[1] * dfDeltaX;
//      break;
    case GDT_UInt16:
        BILINEAR_INTERPOLATION(unsigned short);
        break;
    case GDT_Int16:
        BILINEAR_INTERPOLATION(short);
        break;
    case GDT_UInt32:
        BILINEAR_INTERPOLATION(unsigned int);
        break;
    case GDT_Int32:
        BILINEAR_INTERPOLATION( int);
        break;
    case GDT_Float32:
        BILINEAR_INTERPOLATION(float);
        break;
    case GDT_Float64:
        BILINEAR_INTERPOLATION(double);
        break;
    default:
        BILINEAR_INTERPOLATION(BYTE);   
        break;
    }
    
    return NOVALUE;
}
double  CImageBase::GetBandVal(double col,double row,int band)
{
    if (col > GetCols()-1 || row > GetRows()-1 || col<0 || row<0 || band<0 || band > GetBands() ) return NOVALUE;

    if (!m_pImgBuf && !malloc_data_buf(m_nBufCols, m_nBufRows, GetBands(),m_pImgSet->GetByte()))     return NOVALUE;

    adjust_data_buf(col, row);

    return GetBandBufVal(col-m_stBufCol,row-m_stBufRow,band);
}

BYTE    CImageBase::GetBandVal8(int col, int row, int band)
{
    static int colorTab[65536]; static int band_last = -1;

    if (col > GetCols() - 1 || row > GetRows() - 1 || col<0 || row<0 || band<0 || band > GetBands()) return 0;
    if (!m_pImgBuf && !malloc_data_buf(m_nBufCols, m_nBufRows, GetBands(), m_pImgSet->GetByte()))    return 0;
    adjust_data_buf(col, row);

    int datasize = m_pImgSet->GetByte();
    BYTE* pBuf = m_pImgBuf + (band*m_nBufRows*m_nBufCols + (row-m_stBufRow)*m_nBufCols + col-m_stBufCol)*datasize;
    switch (datasize)
    {
    case 2:
    {
              if (band_last != band) { m_pImgSet->Get16to8Tab(colorTab, band);  band_last = band; }
              return (BYTE)colorTab[*(unsigned int*)pBuf];
    }
    default:
        return *pBuf;
    }
}

void   CImageBase::GetGeoXY(double* x, double* y) const
{
    m_pImgSet->GetGeoXY(x, y);
}
void   CImageBase::GetImgXY(double* x, double* y) const
{
    m_pImgSet->GetImgXY(x, y);
}

bool    CImageBase::GetStatistics(int nBandIdx, double *pdfMin, double *pdfMax, double *pdfMean, double *padfStdDev) const{
    return m_pImgSet->GetStatistics(nBandIdx, pdfMin, pdfMax, pdfMean, padfStdDev);
}

bool    CImageBase::GetNoDataValue(int nBandIdx, double* val)
{
    return m_pImgSet->GetNoDataValue(nBandIdx, val);
}
bool    CImageBase::SetNoDataValue(int nBandIdx, double val){
    return m_pImgSet->SetNoDataValue(nBandIdx, val);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//
bool InvGeoTransform(const double *gt_in, double *gt_out)

{
    double  det, inv_det;
    det = gt_in[1] * gt_in[5] - gt_in[2] * gt_in[4];

    if (fabs(det) < 0.000000000000001)
        return false;
    inv_det = 1.0 / det;

    gt_out[1] = gt_in[5] * inv_det;
    gt_out[4] = -gt_in[4] * inv_det;

    gt_out[2] = -gt_in[2] * inv_det;
    gt_out[5] = gt_in[1] * inv_det;

    gt_out[0] = (gt_in[2] * gt_in[3] - gt_in[0] * gt_in[5]) * inv_det;
    gt_out[3] = (-gt_in[1] * gt_in[3] + gt_in[0] * gt_in[4]) * inv_det;

    return true;
}
