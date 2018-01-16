
#include "LxString.hpp"

#include "OGRBase.h"

#include <ogrsf_frmts.h>

struct ogr_dll_register
{
    ogr_dll_register()
    {
        OGRRegisterAll();
    }
} g_ogr_dll_register;

//////////////////////////////////////////////////////////////////////////////////////////////////
//COGRFeatureBase
#ifndef _VERTICE
#define _VERTICE
typedef struct tagvertice{
    OGRPolygon      *pPolygon;
    OGRLineString   *pLine;
    int             lineIdx;
    OGRPoint        *pPoint;
    int             pointIdx;
}vertice;
#endif
#define VERTICE_PTR(x)  ((vertice*)x)
#define FEATURE_PTR(x)  ((OGRFeature*)x)

COGRFeatureBase::COGRFeatureBase()
{
//  m_pOGRFeature = new OGRFeature[1];
    m_pOGRFeature = NULL;
    m_pOGRFeatureAttach = NULL;
    m_pVertice = new vertice;
}
COGRFeatureBase::~COGRFeatureBase()
{
    if (m_pOGRFeature)  OGRFeature::DestroyFeature((OGRFeature*)m_pOGRFeature); m_pOGRFeature = NULL;
    if (m_pVertice) delete (vertice*)m_pVertice; m_pVertice = NULL;
}

int COGRFeatureBase::GetFieldCount()
{
    return ((OGRFeature*)GetActive())->GetFieldCount();
}

int COGRFeatureBase::GetFieldIndex(const char * pszName)
{
    return ((OGRFeature*)GetActive())->GetFieldIndex(pszName);
}

const char*         COGRFeatureBase::GetFieldName(int i)
{
    OGRFieldDefn *pDefn = ((OGRFeature*)GetActive())->GetFieldDefnRef(i);
    return pDefn->GetNameRef();
}

OGRFieldDataType    COGRFeatureBase::GetFieldType(int idx){
    return (OGRFieldDataType)((OGRFeature*)GetActive())->GetFieldDefnRef(idx)->GetType();
}
int                 COGRFeatureBase::GetFieldAsInteger(int i){
    return ((OGRFeature*)GetActive())->GetFieldAsInteger(i);
}
double              COGRFeatureBase::GetFieldAsDouble(int i){
    return ((OGRFeature*)GetActive())->GetFieldAsDouble(i);
}
const char         *COGRFeatureBase::GetFieldAsString(int i){
    return ((OGRFeature*)GetActive())->GetFieldAsString(i);
}
const int          *COGRFeatureBase::GetFieldAsIntegerList(int i, int *pnCount){
    return ((OGRFeature*)GetActive())->GetFieldAsIntegerList(i,pnCount);
}
const double       *COGRFeatureBase::GetFieldAsDoubleList(int i, int *pnCount){
    return ((OGRFeature*)GetActive())->GetFieldAsDoubleList(i,pnCount);
}
char              **COGRFeatureBase::GetFieldAsStringList(int i){
    return ((OGRFeature*)GetActive())->GetFieldAsStringList(i);
}
BYTE              * COGRFeatureBase::GetFieldAsBinary(int i, int *pnCount){
    return ((OGRFeature*)GetActive())->GetFieldAsBinary(i,pnCount);
}
int                 COGRFeatureBase::GetFieldAsDateTime(int i,
    int *pnYear, int *pnMonth, int *pnDay,
    int *pnHour, int *pnMinute, int *pnSecond,
    int *pnTZFlag){
    return ((OGRFeature*)GetActive())->GetFieldAsDateTime(i,
        pnYear, pnMonth,pnDay,
        pnHour,pnMinute,pnSecond,
        pnTZFlag);
}

int                 COGRFeatureBase::GetFieldAsInteger(const char *pszFName){
    return ((OGRFeature*)GetActive())->GetFieldAsInteger(pszFName);
}
double              COGRFeatureBase::GetFieldAsDouble(const char *pszFName){
    return ((OGRFeature*)GetActive())->GetFieldAsDouble(pszFName);
}
const char         *COGRFeatureBase::GetFieldAsString(const char *pszFName){
    return ((OGRFeature*)GetActive())->GetFieldAsString(pszFName);
}
const int          *COGRFeatureBase::GetFieldAsIntegerList(const char *pszFName, int *pnCount){
    return ((OGRFeature*)GetActive())->GetFieldAsIntegerList(pszFName,pnCount);
}
const double       *COGRFeatureBase::GetFieldAsDoubleList(const char *pszFName, int *pnCount){
    return ((OGRFeature*)GetActive())->GetFieldAsDoubleList(pszFName, pnCount);
}
char              **COGRFeatureBase::GetFieldAsStringList(const char *pszFName){
    return ((OGRFeature*)GetActive())->GetFieldAsStringList(pszFName);
}

OGRwkbType          COGRFeatureBase::GetGeometryType()
{
    OGRFeature* pFeature = (OGRFeature*)GetActive();
    return (OGRwkbType)wkbFlatten((pFeature->GetGeometryRef())->getGeometryType());
}

const char*         COGRFeatureBase::GetGeometryTypeToName(OGRwkbType type)
{
    return OGRGeometryTypeToName((OGRwkbGeometryType)type);
}
const char*         COGRFeatureBase::GetGeometryTypeToName()
{
    return GetGeometryTypeToName(GetGeometryType());
}

const char*         COGRFeatureBase::GetName()
{
    OGRFeatureDefn* pDefn = ((OGRFeature*)GetActive())->GetDefnRef();
    return pDefn->GetName();
}

int                 COGRFeatureBase::GetGeometryVerticesNum()
{   
    switch (GetGeometryType()){
    case wkb_Point:
        return 1;
    case wkb_LineString:
    {
                           OGRLineString*   pGeometry = (OGRLineString*)((OGRFeature*)GetActive())->GetGeometryRef();
                           return pGeometry->getNumPoints();
    }       
    case wkb_Polygon:
    {
                        OGRPolygon* pGeometry = (OGRPolygon*)((OGRFeature*)GetActive())->GetGeometryRef();
                        int num = 0;
                        
                        num += pGeometry->getExteriorRing()->getNumPoints();
                        for (int i = 0; i < pGeometry->getNumInteriorRings(); i++){
                            num += pGeometry->getInteriorRing(i)->getNumPoints();
                        }
                        return num;
    }
    }
    return 0;
}

void            COGRFeatureBase::ResetVerticesReading()
{
    memset(m_pVertice, 0, sizeof(vertice));
    switch (GetGeometryType()){
    case wkb_Point:
    {
                      VERTICE_PTR(m_pVertice)->pPoint = (OGRPoint *)FEATURE_PTR(GetActive())->GetGeometryRef();
                      break;
    }
    case wkb_LineString:
    {
                      VERTICE_PTR(m_pVertice)->pLine = (OGRLineString *)FEATURE_PTR(GetActive())->GetGeometryRef();
                      break;
    }
    case wkb_Polygon:
    {
                      VERTICE_PTR(m_pVertice)->pPolygon = (OGRPolygon *)FEATURE_PTR(GetActive())->GetGeometryRef();
                      VERTICE_PTR(m_pVertice)->pLine = VERTICE_PTR(m_pVertice)->pPolygon->getExteriorRing();
                      break;
    }
    default:
        break;
    }
}

inline void     GetPoint(OGRPoint *pPoint, double* x, double* y, double* z){
    if (x) *x = pPoint->getX();
    if (y) *y = pPoint->getY();
    if (z) *z = pPoint->getZ();
}

inline bool GetCurPointPtr(OGRPoint* pt,int idx, OGRLineString* pLine){
    if (idx >= pLine->getNumPoints()) return false;
    ;
    pLine->getPoint(idx, pt);
    return true;
}

inline OGRLineString* GetCurLine(int idx, OGRPolygon* pPolygon){
    if (idx >= pPolygon->getNumInteriorRings() + 1) return NULL;
    if (idx == 0) return pPolygon->getExteriorRing();
    return pPolygon->getInteriorRing(idx-1);
}

bool            COGRFeatureBase::GetNextVertices(double* x, double* y, double* z)
{
    OGRPoint pt;
    if (VERTICE_PTR(m_pVertice)->pPolygon){
        OGRLineString* pLine;
        
        if (!GetCurPointPtr(&pt,VERTICE_PTR(m_pVertice)->pointIdx, VERTICE_PTR(m_pVertice)->pLine)) {
            VERTICE_PTR(m_pVertice)->lineIdx++;
            pLine = GetCurLine(VERTICE_PTR(m_pVertice)->lineIdx, VERTICE_PTR(m_pVertice)->pPolygon);
            if (!pLine) return false;           
            VERTICE_PTR(m_pVertice)->pointIdx = 0;
            VERTICE_PTR(m_pVertice)->pLine = pLine;
            pLine->getPoint(0, &pt);
        }
        VERTICE_PTR(m_pVertice)->pPoint = &pt;
    }
    else if (VERTICE_PTR(m_pVertice)->pLine){       
        
        if (!GetCurPointPtr(&pt,VERTICE_PTR(m_pVertice)->pointIdx, VERTICE_PTR(m_pVertice)->pLine)) return false;
        VERTICE_PTR(m_pVertice)->pPoint = &pt;
    }
    else if (VERTICE_PTR(m_pVertice)->pPoint){
        if (VERTICE_PTR(m_pVertice)->pointIdx >= 1) return false;
    }
    else
    {
        return false;
    }

    GetPoint(VERTICE_PTR(m_pVertice)->pPoint, x, y, z);
    VERTICE_PTR(m_pVertice)->pointIdx++;
    return true;
}

void COGRFeatureBase::Attach(void* pFeature){
    m_pOGRFeatureAttach = pFeature;
    ResetVerticesReading();
}

// bool COGRFeatureBase::BeginCreate(OGRwkbType eType){
//  switch (eType)
//  {
//  default:
//      break;
//  }
// }
// 
// bool COGRFeatureBase::AddPoint(int num,double* x,double* y,double* z /* = NULL */)
// {
//  ((OGRFeature*)GetActive())->
// }

//////////////////////////////////////////////////////////////////////////////////////////////////
//COGRLayerBase
//create in mem firstly
COGRLayerBase::COGRLayerBase(){
//  m_pOGRLayer = new OGRLayer[1];
    m_pOGRLayer = NULL;
    m_pOGRLayerAttach = NULL;
}
COGRLayerBase::~COGRLayerBase(){
    if (m_pOGRLayer)  delete (OGRLayer*)m_pOGRLayer; m_pOGRLayer = NULL;
}

void COGRLayerBase::SetSpatialFilterRect(double dfMinX, double dfMinY, double dfMaxX, double dfMaxY)
{
    ((OGRLayer*)GetActive())->SetSpatialFilterRect(dfMinX, dfMinY, dfMaxX, dfMaxY);
}

bool COGRLayerBase::SetAttributeFilter(const char * filter)
{
    if (((OGRLayer*)GetActive())->SetAttributeFilter(filter) == CE_Failure){
        return false;
    }
    return true;
}

void COGRLayerBase::ResetReading()
{
    ((OGRLayer*)GetActive())->ResetReading();
}

bool COGRLayerBase::GetNextFeature(COGRFeatureBase& feature)
{
    OGRFeature * pFeature = ((OGRLayer*)GetActive())->GetNextFeature();
    if (!pFeature)  return false;
    feature.Attach(pFeature);
    return true;
}

int  COGRLayerBase::GetFeatureCount(bool bForce )
{
    return ((OGRLayer*)GetActive())->GetFeatureCount(bForce ? 1 : 0);
}

const char* COGRLayerBase::GetName()
{
    return ((OGRLayer*)GetActive())->GetName();
}

bool        COGRLayerBase::CreateField(const char* lpstrFieldName, OGRFieldDataType type)
{
    OGRFieldDefn oField(lpstrFieldName, (OGRFieldType)type);
    if (((OGRLayer*)GetActive())->CreateField(&oField) != OGRERR_NONE) return false;
    return true;    
}

void        COGRLayerBase::Create(void* pLayer){
    if (m_pOGRLayer)  delete (OGRLayer*)m_pOGRLayer; m_pOGRLayer = NULL;
    m_pOGRLayer = pLayer;
}
//////////////////////////////////////////////////////////////////////////////////////////////////
//COGRDataSetBase
COGRDataSetBase::COGRDataSetBase()
{
    m_pOGRSet = NULL;
}


COGRDataSetBase::~COGRDataSetBase()
{
    if (m_pOGRSet) OGRDataSource::DestroyDataSource((OGRDataSource *)m_pOGRSet); m_pOGRSet = NULL;
}

bool COGRDataSetBase::Open(const char* lpstrPath){
    Reset();

    m_pOGRSet = OGRSFDriverRegistrar::Open(lpstrPath, FALSE);
    if (!m_pOGRSet) {
        LogPrint(ERR_FILE_FLAG, GetLastErrorMsg()); return false;
    }
    SetFilePath(lpstrPath);
    return true;
}

bool COGRDataSetBase::Create(const char* lpstrPath){
    Reset();
    const char* pExt = strrchr(lpstrPath, '.'); pExt = GetDescription(pExt);
    if (!pExt) { LogPrint(ERR_FILE_FLAG, "unsupported file type %s", lpstrPath); return false; }
    OGRSFDriver* poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(pExt);
    m_pOGRSet = poDriver->CreateDataSource(lpstrPath, NULL);
    if (!m_pOGRSet) return false;
    SetFilePath(lpstrPath);
    return true;
}

bool COGRDataSetBase::CreateLayer(COGRLayerBase& layer, const char* lpstrLayerName, OGRwkbType wkbType, const char* lpstrWkt /* = NULL */)
{
    OGRSpatialReference srcSpaRef;  bool bFlag = true;
    if (!lpstrWkt || srcSpaRef.importFromWkt((char **)&lpstrWkt) == CE_Failure)
        bFlag = false;
    OGRLayer*  player = ((OGRDataSource *)m_pOGRSet)->CreateLayer(lpstrLayerName, bFlag ? &srcSpaRef : NULL, (OGRwkbGeometryType)wkbType);
    if (!player) return false;
    layer.Create(player);
    return true;
}

void COGRDataSetBase::Close(){
    Reset();
}
int  COGRDataSetBase::GetLayerCount(){
    return ((OGRDataSource *)m_pOGRSet)->GetLayerCount();
}
bool COGRDataSetBase::InitLayerByIdx(COGRLayerBase& layer, int idx)
{
    OGRLayer*  player = ((OGRDataSource *)m_pOGRSet)->GetLayer(idx);
    if (!player) return false;
    layer.Attach(player);
    return true;
}
bool COGRDataSetBase::InitLayerByName(COGRLayerBase& layer, const char* layerName)
{
    OGRLayer*  player = ((OGRDataSource *)m_pOGRSet)->GetLayerByName(layerName);
    if (!player) return false;
    layer.Attach(player);
    return true;
}

const char* COGRDataSetBase::GetLastErrorMsg()  { 
    return CPLGetLastErrorMsg(); 
}

void    COGRDataSetBase::Reset(){
    memset(m_strOGRPath, 0, 512);
    if (m_pOGRSet) OGRDataSource::DestroyDataSource((OGRDataSource *)m_pOGRSet);    m_pOGRSet = NULL;
}

const char* COGRDataSetBase::GetDescription(const char* lpstrExt)
{
    static const char* strExt[] = { "shp" };
    static const char* strDescrip[] = { "ESRI Shapefile" };

    char suffix[10];    if (*lpstrExt == '.') strcpy(suffix, lpstrExt + 1); else strcpy(suffix, lpstrExt);
    strlwr(suffix);

    for (int i = 0; i < sizeof(strExt) / sizeof(const char*); i++){
        if (!strcmp(strExt[i], suffix)) {
            return strDescrip[i];
        }
    }
    return NULL;
}
