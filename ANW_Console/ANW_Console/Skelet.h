#include "memory.h"

unsigned unsigned char mask=0x80;


int skelet (unsigned unsigned char *b, int ns, int ls );
int readpix (unsigned char *b, unsigned char *bx, int ls, int i, int j );
void outbit (unsigned char *b, int ls, int i, int j, int zn );
int checkconfig (unsigned char *b, unsigned char *bx, int ns, int ls, int i, int j );
void readconfig (unsigned char *b, unsigned char *bx, int ns, int ls, int i, int j,
				  int mconf[] );


int skelet (unsigned char *b, int ns, int ls )
  {
	unsigned char *bx;
	int i,j,js,nj,jn,in;
	int prism,prost;
	int pix,pixneib;

	  bx = (unsigned char*)malloc (ns*ls);
	  if(bx==NULL)
			 return 1;
	  memset(bx,'\0',ns*ls);
	  nj=8*ls;
	  prism=1;
	  while(prism)
		{
		  prism=0;
		  for(js=0;js<8;js+=2)
			{
			  for(i=0;i<ns;i++)
				for(j=0;j<nj;j++)
				  {
					pix=readpix(b,bx,ls,i,j);
					if(pix==1)
					  {
						switch(js)
						  {
							case 0 : in=i; jn=j+1;
									 break;
							case 2 : in=i-1; jn=j;
									 break;
							case 4 : in=i; jn=j-1;
									 break;
							case 6 : in=i+1; jn=j;
									 break;
						   }
						if(in<0 || in==ns || jn<0 || jn==nj)
						  pixneib=0;
						else
						  pixneib=readpix(b,bx,ls,in,jn);
						if(pixneib==0)
						  {
							prost=checkconfig(b,bx,ns,ls,i,j);
							outbit(bx,ls,i,j,1);
							if(prost);
							else
							  {
								outbit(b,ls,i,j,0);
							   }
							prism=1;
						   }
					   } 
				   } 
			  if(prism)
				for(i=0;i<ns;i++)
				  for(j=0;j<nj;j++)
					{
					  pix=readpix(b,bx,ls,i,j);
					  if(pix==3)
						outbit(bx,ls,i,j,0);
					 }
			 } 
		 }
	  free(bx);
	  return 0;
   }  

int readpix (unsigned char *b,unsigned char *bx, int ls, int i, int j )
  {
	unsigned char c,cx,m;
	int jb,jm,pix;
	unsigned off;
	  jb=j/8;
	  jm=j%8;
	  m=mask>>jm;
	  off=i*ls+jb;
	  c=*(b+off);
	  cx=*(bx+off);
	  c=c&m;
	  cx=cx&m;
	  if(c&cx)
		pix=2;
	  else
		if((~c)&cx)
		  pix=3;
		else
		  if(c&(~cx))
			pix=1;
		  else
			pix=0;
	  return pix;
   }    

void outbit (unsigned char *b, int ls, int i, int j, int zn )
  {
	unsigned char c,m;
	int jb,jm;
	unsigned off;
	  jb=j/8;
	  jm=j%8;
	  off=i*ls+jb;
	  c=*(b+off);
	  m=mask>>jm;
	  if(zn) c=c | m;
		else c=c & (~m);
	  *(b+off)=c;
   }


 char *stconf[]={"0AAA0BBB","AA0BBB0A","0AAAAA02",
				"AAAA020A","AA020AAA","020AAAAA"};

int checkconfig (unsigned char *b, unsigned char *bx, int ns, int ls, int i, int j )
  {
	int mconf[8];
	int istconf,jj;
	int prA,prB,pr2;
	int z,prost,przap;
	unsigned char c;
	  readconfig(b,bx,ns,ls,i,j,mconf);
	  prost=0;
	  for(istconf=0;istconf<6;istconf++)
		{
		  prA=0;
		  prB=0;
		  pr2=0;
		  przap=0;
		  for(jj=0;jj<8;jj++)
			{
			  c=*(stconf[istconf]+jj);
			  z=mconf[jj];
			  if(c=='A' && z>0) prA=1;
			  if(c=='B' && z>0) prB=1;
			  if(c=='2' && (z==1 || z==2)) pr2=1;
			  if(c=='0' && z!=0) przap=1;
			 }
		  if((prA && prB || prA && pr2) && przap==0)
			{
			  prost=1;
			  break;
			 }
		 }
	  return prost;
   }  

void readconfig (unsigned char *b,unsigned char *bx, int ns, int ls, int i, int j,
				  int mconf[] )
  {
	int js,in,jn;
	int pixn;
	int nj;
	  nj=8*ls;
	  for(js=0;js<8;js++)
		{
		  switch(js)
			{
			  case 0 : in=i; jn=j+1;
					   break;
			  case 1 : in=i-1; jn=j+1;
					   break;
			  case 2 : in=i-1; jn=j;
					   break;
			  case 3 : in=i-1; jn=j-1;
					   break;
			  case 4 : in=i; jn=j-1;
					   break;
			  case 5 : in=i+1; jn=j-1;
					   break;
			  case 6 : in=i+1; jn=j;
					   break;
			  case 7 : in=i+1; jn=j+1;
					   break;
			 }
		  if(in<0 || in==ns || jn<0 || jn==nj)
			pixn=0;
		  else
			pixn=readpix(b,bx,ls,in,jn);
		  mconf[js]=pixn;
		 }
   } 