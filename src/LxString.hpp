////LxString.hpp
/********************************************************************
	LxString
	created:	2013/10/01
	author:		LX 
	purpose:	This file is for LxString function
*********************************************************************/
#if !defined LxString_hpp__LX_2013_10_1
#define LxString_hpp__LX_2013_10_1

#include "PrjLog.hpp"

//-------------size: bit of *p
inline void	BinaryString(char* str,const void* p,size_t size)
{
	int j;
	const char* pS = (const char*) p;
	size_t cnt=0;

	while (1){
		char a = 0x01;

		for (j=0;j<8;j++){	//sizeof(char)
			*str++ = *pS&a? '1':'0';
			a<<=1;

			cnt++;
			if( cnt>=size ) return;
		}
		pS++;
	}
			
}

inline void BinaryString2int(int* p,LPCSTR str,size_t size)
{
	*p=0;
	for (int i=size-1; i>=0; i--){
		*p <<= 1;
		if( str[i] != '0' ) *p += 1;
	}

}

inline short strchrPos (
					   const char * string,
					   int ch
					   )
{
	short pos = 1;
	while (*string && *string != (char)ch)
	{ string++; pos++; }
	
	if (*string == (char)ch)
		return(pos);
	return(0);
}

inline short strrchrPos (
						const char * string,
						int ch
						)
{
	short pos = 0;
	char *start = (char *)string;
	
	while (*string++)                       /* find end of string */
		pos++;
	/* search towards front */
	while (--string != start && *string != (char)ch)
		pos--;
	
	if (*string == (char)ch)                /* char found ? */
		return( pos+1 );
	
	return(0);
}

inline short strccpy(
					 const char * src,
					 char *	dst,
					 int ch
					 )
{
	short len = 0;
	char *start = (char *)src;

	while (*src && *src != (char)ch)
	{ src++; len++; }

	if (*src == (char)ch)
	{
		memcpy(dst,start,len);
		dst[len] = 0;
		return(len);
	}
	return(0);
}

inline bool strrcnct(char* dst,const char c,const char* src){
	char* pT = strrchr(dst,c);	
	if(!pT) return false;

	strcpy(pT,src);
	return true;
}

inline char* AddEndSlash( char *strPN ){
	if( strlen(strPN)<1 ) return NULL;
	char* pS = strPN + strlen(strPN)-1;
	if(*pS == '/'||*pS == '\\') return pS;
	else if(strPN[1]!=':' && strPN[1]!='\\'){
		pS++;	*pS = '/'; *(pS+1) = '\0';
	}
	else { pS++;	*pS = '\\'; *(pS+1) = '\0'; }
	return pS;
}

inline bool RelativePath2AbsolutePath(const char* rlt,const char* current_dir,char* abs)
{
	char path[512];	strcpy(path,current_dir);
	int len = strlen(path);	if(len<1) { len = 1; }
	char* pD = path+len-1;
	if( *pD=='/'||*pD=='\\' ) *pD=0;
	
	const char* pS = rlt;	
	if( *pS=='.' && (*(pS+1)=='\\'||*(pS+1)=='/') ) { sprintf(abs,"%s%s",path,pS+1); return true;	}
	
	while( *pS=='.' && *(pS+1)=='.' && (*(pS+2)=='\\'||*(pS+2)=='/')) 
	{
		pD = strrchr(path,'/');		if(pD) pD = strrchr(path,'\\');	if(!pD) { *abs=0; return false; }
		*pD = 0;
		pS += 3;
	}

	*pD = '/';
	sprintf(abs,"%s%s",path,pS);	
	return true;
}

inline bool CheckChar(const char c){
	if(c == '\n') return false;
	if(c == '\r') return false;
	if(c == '\0') return false;
	return true;
}

//return end position of long string
inline int	MatchString(LPCSTR	lpLongString,LPCSTR lpShortString){
	if( strlen(lpLongString) < strlen(lpShortString) ) return 0;
	
	const char *pL,*pS;	int nMark,shortMark;
	pL=lpLongString;	pS=lpShortString;
	nMark=0;	shortMark = 0;

	char headchar = *lpShortString;	int headcnt = 0;
	while( CheckChar(*pS)&&(headchar==*pS++) ) 
		headcnt++;
	
	pS = lpShortString;
	while( CheckChar(*pL) ){
		if( *pL == *pS ) {
			pS++;	shortMark++;				//nMark++;	
		}
		else	{ 
			if( headcnt!=shortMark || *pL!=headchar ) { pS=lpShortString; shortMark=0; }
		}	//nMark = 0;
		pL++;		nMark++;
		if( !CheckChar(*pS) ) 
			return nMark;
	}
	return 0;
}

inline int	PickString(LPCSTR	lpLongString,LPCSTR lpShortString){
	if( strlen(lpLongString) < strlen(lpShortString) ) return 0;
	int nMark;
	char strLong[1024],strShort[512];	bool bAddSpace_L=false;
	if(*lpLongString== ' ') strcpy(strLong,lpLongString);	else { sprintf(strLong," %s",lpLongString); bAddSpace_L=true; }
	if(*lpShortString== ' ') strcpy(strShort,lpShortString);	else { sprintf(strShort," %s",lpShortString); }
	
	nMark = MatchString(strLong,strShort);
	if(bAddSpace_L) nMark -= 1;
	if(nMark<0) nMark =0;
	
	return nMark;
}

inline bool ReplaceString( char* strLine,const char* strName,const char* strReplace )
{
	int nTmp=0;
	nTmp = MatchString(strLine,strName);

	if (nTmp==0) return false;
	int nNamLen = strlen(strName);	nNamLen = nTmp-nNamLen;
	if(nNamLen<0) return false;
	
	char strH[512];	memcpy(strH,strLine,nNamLen*sizeof(char));
	memcpy(strH+nNamLen,strReplace,strlen(strReplace)*sizeof(char));	nNamLen += strlen(strReplace);
	strcpy(strH+nNamLen,strLine+nTmp);
	strcpy(strLine,strH);
	
	return true;
}

inline bool GetPrivateProfilePath(char* filepath,LPCSTR lpstrFileName)
{
	filepath[0] = 0;
	char execname[256]="";  readlink( "/proc/self/exe",execname,256 );
	Dos2Unix(execname);
	char	strFileName[128]="";

	char* pS = NULL;	pS = strrchr(execname,'.');	if(!pS) pS = execname+strlen(execname);	strcpy(pS,".ini");
//	strcpy(filepath,execname);
	pS = strrchr(execname,'/');		if(!pS) return false;
	if( lpstrFileName )		{	strcpy(pS+1,lpstrFileName);	strcpy(strFileName,lpstrFileName);	}
	else strcpy(strFileName,pS+1);

	while(1)
	{
		if (IsExist(execname)) { strcpy(filepath,execname); return true; }
		*pS = 0;
		pS = strrchr(execname,'/');		if(!pS) return false; 
		strcpy(pS+1,strFileName);	 

	}
	
	return false;
}

inline DWORD	GetPrivateProfileStringE( 
										LPCSTR lpAppName, 
										LPCSTR lpKeyName,	
										LPCSTR lpDefault, 
										LPSTR  lpReturnedString,
										DWORD  nSize,
										LPCSTR lpFileName
										)
{
	DWORD nLen;	LPCSTR lpCopy;
	lpCopy = lpDefault;
	
	char strLine[1024];	
	if(IsExist(lpFileName)&&*lpAppName&&*lpKeyName){
		FILE* fp = fopen(lpFileName,"rt");
		int nMark=0,nt;
		bool bTmp=false;
		while(!feof(fp)){
			memset(strLine,0,sizeof(char)*1024);	fgets(strLine,1024,fp);
			if(bTmp)	if( (nMark = PickString(strLine,lpKeyName)) ) { if(strLine[nMark]==' '||strLine[nMark]=='\t'||strLine[nMark]=='=') { break;} }
			if((nMark=PickString(strLine,"[")))	{
				if(bTmp) { nMark=0;break;}
				else if(( nt=PickString(strLine+nMark,lpAppName) ) ) { if(strLine[nt+nMark]==' '||strLine[nt+nMark]==']') bTmp = true; }
			}
		}
		fclose(fp);
		if(bTmp&&(nMark)){
			nLen = strlen(strLine);
			while( strLine[nLen-1]=='\n'||strLine[nLen-1]=='\r'||strLine[nLen-1]=='\t'||strLine[nLen-1]==' ' ) { nLen--; strLine[nLen]=0; }
			while(strLine[nMark] == ' '||strLine[nMark] == '\t'||strLine[nMark] == '=') { nMark++; }
			lpCopy=strLine+nMark;
			if(strlen(lpCopy)==0) lpCopy = lpDefault;
		}
	}
	
	nLen = strlen(lpCopy);	
	if( nLen>=(nSize-1) ) { memcpy(lpReturnedString,lpCopy,nSize-1); nLen = nSize-1; }
	else  memcpy(lpReturnedString,lpCopy,nLen);
	*(lpReturnedString+nLen)=0; 
	
	return nLen+1;
	
}

#ifndef WIN32

#define GetPrivateProfileString	GetPrivateProfileStringE
// inline DWORD	GetPrivateProfileString( 
// 								LPCSTR lpAppName, 
// 								LPCSTR lpKeyName,	
// 								LPCSTR lpDefault, 
// 								LPSTR  lpReturnedString,
// 								DWORD  nSize,
// 								LPCSTR lpFileName
// 								)
// {
// 	DWORD nLen;	LPCSTR lpCopy;
// 	lpCopy = lpDefault;
// 	
// 	char strLine[1024],strVal[512];	strVal[0]=0;
// 	if(IsExist(lpFileName)&&*lpAppName&&*lpKeyName){
// 		FILE* fp = fopen(lpFileName,"rt");
// 		int nMark=0,nt;
// 		bool bTmp=false;
// 		while(!feof(fp)){
// 			fgets(strLine,1024,fp);
// 			if(bTmp)	if( (nMark = PickString(strLine,lpKeyName)) ) { if(strLine[nMark]==' '||strLine[nMark]=='=') { nt=PickString(strLine+nMark,"="); nMark += nt; break;} }
// 			if((nMark=PickString(strLine,"[")))	{
// 				if(bTmp) { nMark=0;break;}
// 				else if(( nt=PickString(strLine+nMark,lpAppName) ) ) { if(strLine[nt+nMark]==' '||strLine[nt+nMark]==']') bTmp = true; }
// 			}
// 		}
// 		fclose(fp);
// 		if(bTmp&&(nMark)){
// 			if(strLine[nMark] == ' '||strLine[nMark] == '=') sscanf(strLine+nMark+1,"%s",strVal); else sscanf(strLine+nMark,"%s",strVal);
// 			lpCopy=strVal;
// 			if(strlen(lpCopy)==0) lpCopy = lpDefault;
// 		}
// 	}
// 	
// 	nLen = strlen(lpCopy);	
// 	if( nLen>=nSize ) { memcpy(lpReturnedString,lpCopy,nSize); nLen = nSize; }
// 	else  memcpy(lpReturnedString,lpCopy,nLen);
// 	*(lpReturnedString+nLen)=0; 
// 	
// 	return nLen;
// 	
// }
// 
#endif

#endif // LxString_hpp__LX_2013_10_1