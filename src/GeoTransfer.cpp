
#include <gdal_priv.h>
#include <ogr_spatialref.h>

#include "LxString.hpp"
#include "GeoTransfer.h"

#define SPGC_PI                 3.1415926535897932384626433832795e0
#define SPGC_R2D                57.295779513082320876798154814105e0
#define SPGC_D2R                0.017453292519943295769236907684886e0

enum ELPS_TYPE
{
    ET_WGS84 = 0,
    ET_BEIJING54 = 1,
    ET_XIAN80 = 2,
    ET_CN2000 = 3,
    ET_CUSTOM = 4,
};

#define SEMIMAJOR_WGS84         6378137.0   
#define SEMIMINOR_WGS84         6356752.31424517929

#define SEMIMAJOR_BEIJING54     6378245.0   
#define SEMIMINOR_BEIJING54     6356863.018773047267

#define SEMIMAJOR_XIAN80        6378140.0   
#define SEMIMINOR_XIAN80        6356755.2881575287

#define SEMIMAJOR_CN2000        6378137.0   
#define SEMIMINOR_CN2000        6356752.3141

const char*     gcs_tag[] = { "WGS 84", "GCS_Beijing_1954", "GCS_Xian_1980", "GCS_China_Geodetic_Coordinate_System_2000","GCS_USER" };//"GCS_Geographic Coordinate System",
const char*     datum_key_tag[] = { "WGS", "BEIJING", "XIAN", "China_2000", "USER" };
const char*     datum_tag[] = { "D_WGS84", "D_BEIJING_1954", "D_XIAN 1980", "China_2000", "USER" };
const char*     ellipsoidName_tag[] = {"WGS84","Krassovsky","Xian_1980","CGCS2000",""};
const double    elps_a[] = { SEMIMAJOR_WGS84, SEMIMAJOR_BEIJING54, SEMIMAJOR_XIAN80, SEMIMAJOR_CN2000 };
const double    elps_b[] = { SEMIMINOR_WGS84, SEMIMINOR_BEIJING54, SEMIMINOR_XIAN80, SEMIMINOR_CN2000 };

enum PROJECTION_TYPE
{
    GEODETIC_LONLAT = 0,   //经纬度
    TM_PROJECTION = 1,   //横轴墨卡托投影   
    GAUSS_PROJECTION = 2,   //高斯－克吕格投影  
    UTM_PROJECTION = 3,   //通用横轴墨卡托投影
    LAMBERT_PROJECTION = 4,   //兰波特等角投影
    GEODETIC_CENTRXYZ = 5,   //地心坐标    
};

const char*     proj_key_tag[] = { "", SRS_PT_TRANSVERSE_MERCATOR, SRS_PT_GAUSSSCHREIBERTMERCATOR,"UTM","Lambert","" };

enum ELEVATION_TYPE
{
    GEODETIC_HEI = 0,
    GEOID_DATA = 1,
};


//////////////////////////////////////////////////////////////////////////
//
inline int      UTM_Zone(float lon){
    return int((lon + 3) / 6 + 0.5) + 30;
}

inline int      UTM_CentralMeridian(int nZone){
    return  6 * (nZone - 30) - 3;
}

inline int      UTM_CentralMeridian(float lon)
{
    return UTM_CentralMeridian(UTM_Zone(lon));
}

inline int      Guass_Zone(float lon, bool bZoneType_6 = true){
    if (bZoneType_6){
        return int((lon + 3) / 6 + 0.5);
    }
    else
    {
        return int(lon  / 3 + 0.5);
    }
}
inline int      Guass_CentralMeridian(int nZone,bool bZoneType_6 = true){
    if (bZoneType_6){
        return 6 * nZone - 3;
    }
    else
    {
        return 3 * nZone;
    }
}

/*
OGRSpatialReference srs1, srs2;
srs1.importFromEPSG(4978);  // 4978 is the EPSG code for geocentric (ECEF)
srs2.SetWellKnownGeogCS("WGS84");
srs2.SetVertCS("NAVD88 height", "North American Vertical Datum 1988");
srs2.SetExtension("VERT_DATUM", "PROJ4_GRIDS", "g2009conus.gtx");

double x, y, z;
x = -2691744.9992481042;    // A point on the ellipsoid
y = -4262401.8609118778;
z = 3894209.4515372305;

// ecef -> geoid
OGRCoordinateTransformation *oct = OGRCreateCoordinateTransformation(&srs1, &srs2);
oct->Transform(1, &x, &y, &z);
printf(" lon,lat,height %lf,%lf,%lf\n", x, y, z);

If it works, it prints the correct value :
lon, lat, height - 122.272778, 37.871667, 32.269611

If the VertCS / Extension is omitted, or if PROJ cannot find the.gtx file, or anything else goes wrong, you simply get the ellipsoidal height(close to 0.0).
*/

struct geocvt_dll_register
{
    geocvt_dll_register()
    {
        char epsg[256] = "";
        char config[256] = "";  readlink("/proc/self/exe", config, 256);
        char* pS = strrpath(config);    if (!pS) pS = config - 1;
        strcpy(pS + 1, "package_gdal.ini");
        GetPrivateProfileString("GeoCvt", "epsg", "data", epsg, 256, config);
        if (!IsExist(epsg)){
            *(pS + 1) = 0;  strcat(config, epsg);
            strcpy(epsg, config);
        }
        
//      char blank[2] = "\0";   char dot[2] = ";";
//      char* env = getenv("GDAL_DATA");
//      char* str;  
//      if (!env) { env = blank; dot[0] = 0; }
//      str = (char*)malloc(strlen(env) + 512);
//      sprintf(str, "GDAL_DATA=%s%s%s", env,dot, epsg);
// 
//      if (putenv(str) == ENOMEM) printf("fail to add env GDAL_DATA !\n");
// 
//      env = getenv("GDAL_DATA");  printf("GDAL_DATA = %s\n", env);
// 
//      free(str);
        if(IsExist(epsg)) CPLSetConfigOption("GDAL_DATA",epsg);
    }
} g_geocvt_dll_register;

class CGeoTransfer
{
public:
    CGeoTransfer()
    {
        m_pCXYZ2LBH = m_pLBH2CXYZ = m_pLBH2Prj = m_pPrj2LBH = NULL;
        m_srcSpaRef.SetWellKnownGeogCS("WGS84");
        UpdateCT();
        m_pszWKT = NULL;
    }
    ~CGeoTransfer()
    {
        ReleaseCT();
        if (m_pszWKT) {
            CPLFree(m_pszWKT); m_pszWKT = NULL;
        }
    }
    void Debug()
    {

    }
    void ReleaseCT()
    {
        ReleaseCT(m_pPrj2LBH);  ReleaseCT(m_pLBH2Prj);
        ReleaseCT(m_pLBH2CXYZ); ReleaseCT(m_pCXYZ2LBH);
    }
    bool ReadGCDFile(const char* pszGcdPath, OGRSpatialReference* pSrcSpaRef){
        FILE* fp = fopen(pszGcdPath, "r");  if (!fp) return false;
        char strLine[1024]; char* pS;
        int elpsType = 0;
        int projType = 0;
        double a, b;
        double Origin_Latitude;
        double Central_Meridian;
        double False_Easting;
        double False_Northing;
        double Scale_Factor;
        int eleType = 0;
        while (!feof(fp))
        {
            if( !fgets(strLine, 1024, fp) ) break;
            pS = strchr(strLine, '=');  if (!pS) continue;
            pS++;
            if (strstr(strLine, "Ellipsoid")){
                elpsType = atoi(pS);
            }
            else if (strstr(strLine, "Projectid")){
                projType = atoi(pS);
            }
            else if (strstr(strLine, "EllipsoA")){
                a = atof(pS);
            }
            else if (strstr(strLine, "EllipsoB")){
                b = atof(pS);
            }
            else if (strstr(strLine, "ProjectScale")){
                Scale_Factor = atof(pS);
            }
            else if (strstr(strLine, "ProjectOriLat")){
                Origin_Latitude = atof(pS);
            }
            else if (strstr(strLine, "ProjectCenMeri")){
                Central_Meridian = atof(pS);
            }
            else if (strstr(strLine, "ProjectFalseNor")){
                False_Northing = atof(pS);
            }
            else if (strstr(strLine, "ProjectFalseEas")){
                False_Easting = atof(pS);
            }
            else if (strstr(strLine, "Elevationid")){
                eleType = atoi(pS);
            }
        }
        fclose(fp);
        SetCvtPar(elpsType, projType, a, b, Origin_Latitude, Central_Meridian, False_Easting, False_Northing, Scale_Factor, eleType, pSrcSpaRef);
        return true;
    }
    bool Import4WKT(const char* pszWkt)
    {
        OGRSpatialReference srcSpaRef;
        if (srcSpaRef.importFromWkt((char **)&pszWkt) == CE_Failure)
            return false;
        UpdateCT(srcSpaRef);
        return true;
    }
    
    bool Import4GCD(const char* pszGcdPath){
        OGRSpatialReference srcSpaRef;
        if (!ReadGCDFile(pszGcdPath, &srcSpaRef)) return false;
        UpdateCT(srcSpaRef);
        return true;
    }
    bool Export2GCD(const char* pszGcdPath)
    {
        int elpsType = 0;
        int projType = 0;
        double a, b;
        double Origin_Latitude;
        double Central_Meridian;
        double False_Easting;
        double False_Northing;
        double Scale_Factor;
        int eleType = 0;
        GetCvtPar(m_srcSpaRef, &elpsType, &projType, &a, &b, &Origin_Latitude, &Central_Meridian, &False_Easting, &False_Northing, &Scale_Factor, &eleType);
        FILE* fp = fopen(pszGcdPath, "w");  if (!fp) return false;
        
        fprintf(fp, "[GeoCvtInf]\n");
        fprintf(fp, "Ellipsoid=%d\n", elpsType);
        fprintf(fp, "Projectid=%d\n", projType);
        fprintf(fp, "EllipsoA=%lf\n", a);
        fprintf(fp, "EllipsoB=%lf\n", b);
        fprintf(fp, "ProjectScale=%lf\n", Scale_Factor);
        fprintf(fp, "ProjectHemi=1\n");
        fprintf(fp, "ProjectOriLat=%lf\n", Origin_Latitude);
        fprintf(fp, "ProjectCenMeri=%lf\n", Central_Meridian);
        fprintf(fp, "ProjectFalseNor=%lf\n", False_Northing);
        fprintf(fp, "ProjectFalseEas=%lf\n", False_Easting);
        fprintf(fp, "ProjectZoneNo=0\n");
        fprintf(fp, "Elevationid=%d\n",eleType);
        fprintf(fp, "ElevationHeiAdj=\n");

        fclose(fp);
        return true;
    }
    const char* GetWKT(){
        if (m_pszWKT) {
            CPLFree(m_pszWKT); m_pszWKT = NULL;
        }
        m_srcSpaRef.exportToWkt(&m_pszWKT);
        return m_pszWKT;
    }

    bool    SetUTM(int nZone, bool bNorth = true){
        if (m_srcSpaRef.SetUTM(nZone, bNorth ? 1 : 0) == CE_Failure) return false;
        UpdateCT_LBH_PRJ(m_srcSpaRef);
        UpdateLPFUN();
        return true;
    }
    static void GetCvtPar(OGRSpatialReference& srcSpaRef,
        int *elpsType,
        int *projType,
        double *a, double *b,
        double *Origin_Latitude, 
        double *Central_Meridian,
        double *False_Easting,
        double *False_Northing,
        double *Scale_Factor,
        int *eleType){
        const char* pInfo = srcSpaRef.GetAttrValue("DATUM");
        unsigned long i;
        for ( i = 0; i < sizeof(datum_key_tag) / sizeof(datum_key_tag[0]); i++){
            if (strstr(pInfo, datum_key_tag[i])) {
                *elpsType = i;  break;
            }
        }
        OGRErr ret;
        *a = srcSpaRef.GetSemiMajor(&ret);  if (ret != OGRERR_NONE) *a = elps_a[*elpsType>ET_CN2000 ? ET_CN2000 : *elpsType];
        *b = srcSpaRef.GetSemiMinor(&ret);  if (ret != OGRERR_NONE) *b = elps_b[*elpsType>ET_CN2000 ? ET_CN2000 : *elpsType];
        *eleType = GEODETIC_HEI;
        *Origin_Latitude = *Central_Meridian = *False_Easting = *False_Northing = *Scale_Factor = 0.0;

        if (srcSpaRef.IsGeocentric()){
            *projType = GEODETIC_CENTRXYZ;          
            return;
        }
        if (srcSpaRef.IsGeographic()){
            *projType = GEODETIC_LONLAT;
            return;
        }
        pInfo = srcSpaRef.GetAttrValue("PROJCS");//PROJECTION
        for ( i = 1; i < sizeof(proj_key_tag) / sizeof(proj_key_tag[0]) -1; i++){
            if (strstr(pInfo, proj_key_tag[i])) {
                *projType = i;  break;
            }
        }
        *Origin_Latitude = srcSpaRef.GetProjParm(SRS_PP_LATITUDE_OF_ORIGIN);
        *Scale_Factor = srcSpaRef.GetProjParm(SRS_PP_SCALE_FACTOR);
        *Central_Meridian = srcSpaRef.GetProjParm(SRS_PP_CENTRAL_MERIDIAN);
        *False_Easting = srcSpaRef.GetProjParm(SRS_PP_FALSE_EASTING);
        *False_Northing = srcSpaRef.GetProjParm(SRS_PP_FALSE_NORTHING); 
        
    }
    static void SetCvtPar(
        int elpsType,
        int projType,
        double a, double b,
        double Origin_Latitude,  
        double Central_Meridian, 
        double False_Easting,
        double False_Northing,
        double Scale_Factor,
        int eleType,
        OGRSpatialReference* pSrcSpaRef)
    {
        pSrcSpaRef->SetGeogCS(gcs_tag[elpsType], datum_tag[elpsType], ellipsoidName_tag[elpsType],
            elps_a[elpsType], elps_a[elpsType] / (elps_a[elpsType] - elps_b[elpsType]));
        if (projType == GEODETIC_CENTRXYZ){
//          pSrcSpaRef->SetGeocCS(ellipsoidName_tag[elpsType]);
            pSrcSpaRef->SetGeocCS("Geocentric");
            return;
        }
        
        switch (projType)
        {
        case GEODETIC_LONLAT: //经纬度
            break;
        case TM_PROJECTION://横轴墨卡托投影
            pSrcSpaRef->SetTM(Origin_Latitude, Central_Meridian, Scale_Factor, False_Easting, False_Northing);
            break;
        case GAUSS_PROJECTION://高斯－克吕格投影
            pSrcSpaRef->SetGaussSchreiberTMercator(Origin_Latitude, Central_Meridian, Scale_Factor, False_Easting, False_Northing);
            break;
        case UTM_PROJECTION://通用横轴墨卡托投影
            pSrcSpaRef->SetUTM(UTM_Zone(Central_Meridian));
            break;
        case LAMBERT_PROJECTION://兰波特等角投影
            pSrcSpaRef->SetLCC1SP(Origin_Latitude, Central_Meridian, Scale_Factor, False_Easting, False_Northing);
            break;
        default:
            break;
        }
        
        
//      pSrcSpaRef->SetProjParm(SRS_PP_SCALE_FACTOR, Scale_Factor);
//      //////////////////中央经线/////////////////////
//      pSrcSpaRef->SetProjParm(SRS_PP_CENTRAL_MERIDIAN, Central_Meridian / SPGC_D2R);
//      ///////////////东偏///////////////////////
//      pSrcSpaRef->SetProjParm(SRS_PP_FALSE_EASTING, False_Easting);
//      ///////////////北偏///////////////////////
//      pSrcSpaRef->SetProjParm(SRS_PP_FALSE_NORTHING, False_Northing);
//      /////////////////Origin_Latitude/////////////////////
//      pSrcSpaRef->GetProjParm(SRS_PP_LATITUDE_OF_ORIGIN, Origin_Latitude);
    }
public:
    bool    DoNothing(int nCount,
        double *x, double *y, double *z ,
        int *pabSuccess 
        )
    {
        return true;
    }
    bool    Cvt2LBH(int nCount,
        double *x, double *y, double *z,
        int *pabSuccess
        )
    {
        return (this->*m_lpfunCvt2LBH)(nCount, x, y, z, pabSuccess);
    }
    bool    Prj2LBH(int nCount,
        double *x, double *y, double *z ,
        int *pabSuccess )
    {
        return TransformEx(m_pPrj2LBH, nCount, x, y, z, pabSuccess);
    }
    bool    LBH2Prj(int nCount,
        double *x, double *y, double *z ,
        int *pabSuccess )
    {
        return TransformEx(m_pLBH2Prj, nCount, x, y, z, pabSuccess);
    }
    bool    LBH2CXYZ(int nCount,
        double *x, double *y, double *z ,
        int *pabSuccess )
    {
        return TransformEx(m_pLBH2CXYZ, nCount, x, y, z, pabSuccess);
    }
    bool    CXYZ2LBH(int nCount,
        double *x, double *y, double *z ,
        int *pabSuccess )
    {
        return TransformEx(m_pCXYZ2LBH, nCount, x, y, z, pabSuccess);
    }

    bool        IsGeographic() const{
        return m_srcSpaRef.IsGeographic() == 0 ? false : true;
    }
    bool        IsProjected() const{
        return m_srcSpaRef.IsProjected() == 0 ? false : true;
    }
    bool        IsGeocentric() const{
        return m_srcSpaRef.IsGeocentric() == 0 ? false : true;
    }

    const char* GetLastErrorMsg() { return CPLGetLastErrorMsg(); }

    static bool TransformEx(OGRCoordinateTransformation* pCT,
        int nCount,
        double *x, double *y, double *z = NULL,
        int *pabSuccess = NULL){
        if (pCT && pCT->TransformEx(nCount, x, y, z, pabSuccess) != CE_Failure) return true;
        return false;
    }
    static void ReleaseCT(OGRCoordinateTransformation*& pCT){
        if (pCT)
            OCTDestroyCoordinateTransformation((OGRCoordinateTransformationH)pCT);
        pCT = NULL;
    }
protected:
    
    void UpdateLPFUN()
    {
        if (m_srcSpaRef.IsGeographic()){
            m_lpfunCvt2LBH = &CGeoTransfer::DoNothing;
        }
        else if (m_srcSpaRef.IsProjected()){
            m_lpfunCvt2LBH = &CGeoTransfer::Prj2LBH;
        }
        else if (m_srcSpaRef.IsGeocentric()){
            m_lpfunCvt2LBH = &CGeoTransfer::CXYZ2LBH;
        }       
    }
    void UpdateCT(){
        UpdateCT_LBH_CXYZ(m_srcSpaRef);
        UpdateCT_LBH_PRJ(m_srcSpaRef);
        UpdateLPFUN();
    }
    void UpdateCT(OGRSpatialReference& oSpaRef){
        if (m_srcSpaRef.IsSame(&oSpaRef)) return;
        if (m_srcSpaRef.IsSameGeogCS(&oSpaRef)){
            UpdateCT_LBH_PRJ(oSpaRef);
        }
        else{
            UpdateCT_LBH_CXYZ(oSpaRef);
            UpdateCT_LBH_PRJ(oSpaRef);
        }
        
        char* pszWKT = NULL;
        oSpaRef.exportToWkt(&pszWKT);
        char* pszWKTCopy = pszWKT;
        m_srcSpaRef.importFromWkt(&pszWKTCopy);
        CPLFree(pszWKT);
        UpdateLPFUN();
    }

    void UpdateCT_LBH_CXYZ(OGRSpatialReference& srcSpaRef)
    {
        ReleaseCT(m_pLBH2CXYZ);     ReleaseCT(m_pCXYZ2LBH);
        OGRSpatialReference oCXYZSpaRef, oBLHSpaRef;

//      oCXYZSpaRef.importFromEPSG(4978);//WGS84 Geocentric
        oCXYZSpaRef.CopyGeogCSFrom(&srcSpaRef);
        oCXYZSpaRef.SetGeocCS("Geocentric");
        oBLHSpaRef.CopyGeogCSFrom(&srcSpaRef);

        m_pLBH2CXYZ = OGRCreateCoordinateTransformation(&oBLHSpaRef, &oCXYZSpaRef);
        m_pCXYZ2LBH = OGRCreateCoordinateTransformation(&oCXYZSpaRef, &oBLHSpaRef);
    }
    bool UpdateCT_LBH_PRJ(OGRSpatialReference& srcSpaRef){
        ReleaseCT(m_pPrj2LBH);      ReleaseCT(m_pLBH2Prj);
        if (!srcSpaRef.IsProjected()) return false;
        OGRSpatialReference  oBLHSpaRef;
        oBLHSpaRef.CopyGeogCSFrom(&srcSpaRef);
        m_pPrj2LBH = OGRCreateCoordinateTransformation(&srcSpaRef, &oBLHSpaRef);
        m_pLBH2Prj = OGRCreateCoordinateTransformation(&oBLHSpaRef, &srcSpaRef);
        return true;
    }

public:
    friend bool IsSame_SR(const CGeoTransfer* src,const CGeoTransfer* dst){
        return src->m_srcSpaRef.IsSame(&(dst->m_srcSpaRef)) == 0? false : true ;
    }
    friend bool IsSameGeogCS_SR(const CGeoTransfer* src,const CGeoTransfer* dst){
        return src->m_srcSpaRef.IsSameGeogCS(&(dst->m_srcSpaRef)) == 0 ? false : true;
    }
    friend bool IsSameVertCS_SR(const CGeoTransfer* src,const CGeoTransfer* dst){
        return src->m_srcSpaRef.IsSameVertCS(&(dst->m_srcSpaRef)) == 0 ? false : true;
    }
    friend class CCoordinationTransfer;
protected:
    typedef bool (CGeoTransfer::*LPFUNCvt2LBH)(int, double *, double *, double *, int *);
    LPFUNCvt2LBH    m_lpfunCvt2LBH;
private:
    OGRSpatialReference         m_srcSpaRef;
    char*                       m_pszWKT;

//  OGRCoordinateTransformation* m_pCT;
    
    OGRCoordinateTransformation* m_pPrj2LBH;
    OGRCoordinateTransformation* m_pLBH2Prj;
    OGRCoordinateTransformation* m_pLBH2CXYZ;
    OGRCoordinateTransformation* m_pCXYZ2LBH;
};

class CCoordinationTransfer{
public:
    CCoordinationTransfer(){
        m_pDst2Src = m_pSrc2Dst = NULL;
    }
    ~CCoordinationTransfer(){
        ReleaseCT();
    }
    void ReleaseCT(){
        CGeoTransfer::ReleaseCT(m_pSrc2Dst);
        CGeoTransfer::ReleaseCT(m_pDst2Src);
    }
    bool init(const OGRSpatialReference* srcSpaRef, const OGRSpatialReference* dstSpaRef){
        ReleaseCT();
        m_pSrc2Dst = OGRCreateCoordinateTransformation((OGRSpatialReference*)srcSpaRef, (OGRSpatialReference*)dstSpaRef);
        m_pDst2Src = OGRCreateCoordinateTransformation((OGRSpatialReference*)dstSpaRef, (OGRSpatialReference*)srcSpaRef);
        if (!m_pSrc2Dst || !m_pDst2Src) return false;
        return true;
    }
    bool init(const CGeoTransfer* src,const CGeoTransfer* dst){
        return init(&src->m_srcSpaRef, &dst->m_srcSpaRef);
    }
    bool    Src2Dst(int nCount,
        double *x, double *y, double *z,
        int *pabSuccess)
    {
        return CGeoTransfer::TransformEx(m_pSrc2Dst, nCount, x, y, z, pabSuccess);
    }
    bool    Dst2Src(int nCount,
        double *x, double *y, double *z,
        int *pabSuccess)
    {
        return CGeoTransfer::TransformEx(m_pDst2Src, nCount, x, y, z, pabSuccess);
    }
protected:
    OGRCoordinateTransformation* m_pSrc2Dst;
    OGRCoordinateTransformation* m_pDst2Src;
};

////////////////////////////////////////////////////////////////////////////////////////
//CGeoTransfer
CGeoTransferBase::CGeoTransferBase()
{
    m_pGeoTransfer = new CGeoTransfer;
}


CGeoTransferBase::~CGeoTransferBase()
{
    if (m_pGeoTransfer) delete m_pGeoTransfer;  m_pGeoTransfer = NULL;
}

bool CGeoTransferBase::Export2GCD(const char* pszGcdPath) const
{
    return m_pGeoTransfer->Export2GCD(pszGcdPath);
}

bool CGeoTransferBase::Import4WKT(const char* pszWkt)
{
    return m_pGeoTransfer->Import4WKT(pszWkt);
}

bool CGeoTransferBase::Import4GCD(const char* pszGcdPath)
{
    return m_pGeoTransfer->Import4GCD(pszGcdPath);
}
const char* CGeoTransferBase::GetWKT() const
{
    return m_pGeoTransfer->GetWKT();
}

bool CGeoTransferBase::Prj2LBH(int nCount, double *x, double *y, double *z /* = NULL */, int *pabSuccess /* = NULL */)
{
    return m_pGeoTransfer->Prj2LBH(nCount, x, y, z, pabSuccess);
}

bool CGeoTransferBase::LBH2Prj(int nCount, double *x, double *y, double *z /* = NULL */, int *pabSuccess /* = NULL */)
{
    return m_pGeoTransfer->LBH2Prj(nCount, x, y, z, pabSuccess);
}

bool CGeoTransferBase::LBH2CXYZ(int nCount, double *x, double *y, double *z /* = NULL */, int *pabSuccess /* = NULL */)
{
    return m_pGeoTransfer->LBH2CXYZ(nCount, x, y, z, pabSuccess);
}

bool CGeoTransferBase::CXYZ2LBH(int nCount, double *x, double *y, double *z /* = NULL */, int *pabSuccess /* = NULL */)
{
    return m_pGeoTransfer->CXYZ2LBH(nCount, x, y, z, pabSuccess);
}

bool CGeoTransferBase::Cvt2LBH(int nCount, double *x, double *y, double *z /* = NULL */, int *pabSuccess /* = NULL */)
{
    return m_pGeoTransfer->Cvt2LBH(nCount, x, y, z, pabSuccess);
}

bool CGeoTransferBase::LBH2Prj(float* x, float* y, float* z /* = NULL */)
{
    double dx, dy, dz;
    dx = *x;    dy = *y;    if (z) dz = *z;
    if (!m_pGeoTransfer->LBH2Prj(1,&dx,&dy,z?&dz:NULL,NULL)){
        return false;
    }
    *x = (float)dx; *y = (float)dy; if (z) *z = (float)dz;
    return true;
}

bool CGeoTransferBase::Prj2LBH(float* x, float* y, float* z /* = NULL */)
{
    double dx, dy, dz;
    dx = *x;    dy = *y;    if (z) dz = *z;
    if (!m_pGeoTransfer->Prj2LBH(1, &dx, &dy, z ? &dz : NULL,NULL)){
        return false;
    }
    *x = (float)dx; *y = (float)dy; if (z) *z = (float)dz;
    return true;
}

bool CGeoTransferBase::IsGeographic() const{
    return m_pGeoTransfer->IsGeographic();
}

bool CGeoTransferBase::IsProjected() const{
    return m_pGeoTransfer->IsProjected();
}

bool CGeoTransferBase::IsGeocentric() const{
    return m_pGeoTransfer->IsGeocentric();
}

bool CGeoTransferBase::IsSameSpatialReference(CGeoTransferBase* sr) const
{
    return IsSame_SR(m_pGeoTransfer, sr->GetInteriorPtr());
}

bool CGeoTransferBase::IsSameGeogCS(CGeoTransferBase* sr) const
{
    return IsSameGeogCS_SR(m_pGeoTransfer, sr->GetInteriorPtr());
}

bool CGeoTransferBase::IsSameVertCS(CGeoTransferBase* sr) const
{
    return IsSameVertCS_SR(m_pGeoTransfer, sr->GetInteriorPtr());
}

bool CGeoTransferBase::SetUTM_Zone(int nZone, bool bNorth /* = true */)
{
    return m_pGeoTransfer->SetUTM(nZone, bNorth);
}

bool CGeoTransferBase::SetUTM_CenMeri(float lon, bool bNorth /* = true */)
{
    return m_pGeoTransfer->SetUTM(UTM_Zone(lon), bNorth);
}

CCoordinationTransferBase::CCoordinationTransferBase(){
    m_pCoordTransfer = new CCoordinationTransfer;
}
CCoordinationTransferBase::~CCoordinationTransferBase(){
    if (m_pCoordTransfer) delete m_pCoordTransfer; m_pCoordTransfer = NULL;
}

bool CCoordinationTransferBase::init(CGeoTransferBase& src, CGeoTransferBase& dst)
{
    return m_pCoordTransfer->init(src.GetInteriorPtr(), dst.GetInteriorPtr());
}

bool CCoordinationTransferBase::Src2Dst(int nCount, double *x, double *y, double *z, int *pabSuccess)
{
    return m_pCoordTransfer->Src2Dst(nCount, x, y, z, pabSuccess);
}

bool CCoordinationTransferBase::Dst2Src(int nCount, double *x, double *y, double *z /* = NULL */, int *pabSuccess /* = NULL */)
{
    return m_pCoordTransfer->Dst2Src(nCount, x, y, z, pabSuccess);
}
