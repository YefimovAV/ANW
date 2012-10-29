unsigned char *NoiseEdit(unsigned char *k, int wI, int hI);

unsigned char *NoiseEdit(unsigned char *k, int wI, int hI)
{
	int wp = 0, bp = 0;
	for (int i = 0; i < wI*hI; i++)
				{
					if (i == 0)
					{
						wp = 0; bp = 0;
						if(k[i+1] == 1)
							wp++;
						else bp++;
						if(k[i+wI] == 1)
							wp++;
						else bp++;
						if(k[i+1+wI] == 1)
							wp++;
						else bp++;
						if(wp > bp) k[i] = 1;
						if(wp < bp)	k[i] = 0;
					}
					if (i != 0 && i < wI - 1)
					{
						wp = 0; bp = 0;
						if(k[i-1] == 1)
							wp++;
						else bp++;
						if(k[i+1] == 1)
							wp++;
						else bp++;
						if(k[i-1+wI] == 1)
							wp++;
						else bp++;
						if(k[i+wI] == 1)
							wp++;
						else bp++;
						if(k[i+1+wI] == 1)
							wp++;
						else bp++;
						if(wp > bp) k[i] = 1;
						if(wp < bp)	k[i] = 0;				
					}
					if (i == wI-1)
					{
						wp = 0; bp = 0;
						if(k[i-1] == 1)
							wp++;
						else bp++;
						if(k[i-1+wI] == 1)
							wp++;
						else bp++;
						if(k[i+wI] == 1)
							wp++;
						else bp++;
						if(wp > bp) k[i] = 1;
						if(wp < bp)	k[i] = 0;
					}
					if (i != 0 && (i % wI) == 0 && i != wI*hI-wI)
					{
						wp = 0; bp = 0;
						if(k[i-wI] == 1)
							wp++;
						else bp++;
						if(k[i+1] == 1)
							wp++;
						else bp++;
						if(k[i-1+wI] == 1)
							wp++;
						else bp++;
						if(k[i+wI] == 1)
							wp++;
						else bp++;
						if(k[i+1+wI] == 1)
							wp++;
						else bp++;
						if(wp > bp) k[i] = 1;
						if(wp < bp)	k[i] = 0;
					}
					if (i > wI && (i % wI) != 0 && ((i+1) % wI) != 0 && i < wI*hI-wI)
					{
						wp = 0; bp = 0;
						if(k[i-1] == 1)
							wp++;
						else bp++;
						if(k[i+1] == 1)
							wp++;
						else bp++;
						if(k[i-1+wI] == 1)
							wp++;
						else bp++;
						if(k[i+wI] == 1)
							wp++;
						else bp++;
						if(k[i+1+wI] == 1)
							wp++;
						else bp++;
						if(k[i-1-wI] == 1)
							wp++;
						else bp++;
						if(k[i-wI] == 1)
							wp++;
						else bp++;
						if(k[i-wI+1] == 1)
							wp++;
						else bp++;
						if(wp > bp) k[i] = 1;
						if(wp < bp)	k[i] = 0;
					}
					if (i > wI && ((i+1) % wI) == 0 && i != wI*hI-1)
					{
						wp = 0; bp = 0;
						if(k[i-wI] == 1)
							wp++;
						else bp++;
						if(k[i-1] == 1)
							wp++;
						else bp++;
						if(k[i-1-wI] == 1)
							wp++;
						else bp++;
						if(k[i+wI] == 1)
							wp++;
						else bp++;
						if(k[i-1+wI] == 1)
							wp++;
						else bp++;
						if(wp > bp) k[i] = 1;
						if(wp < bp)	k[i] = 0;
					}
					if (i == wI*hI-wI)
					{
						wp = 0; bp = 0;
						if(k[i-wI] == 1)
							wp++;
						else bp++;
						if(k[i+1] == 1)
							wp++;
						else bp++;
						if(k[i+1-wI] == 1)
							wp++;
						else bp++;
						if(wp > bp) k[i] = 1;
						if(wp < bp)	k[i] = 0;
					}
					if (i > wI*hI-wI && i < wI*hI-1)
					{
						wp = 0; bp = 0;
						if(k[i-wI-1] == 1)
							wp++;
						else bp++;
						if(k[i-wI] == 1)
							wp++;
						else bp++;
						if(k[i-1] == 1)
							wp++;
						else bp++;
						if(k[i+1-wI] == 1)
							wp++;
						else bp++;
						if(k[i+1] == 1)
							wp++;
						else bp++;
						if(wp > bp) k[i] = 1;
						if(wp < bp)	k[i] = 0;
					}
					if (i == wI*hI-1)
					{
						wp = 0; bp = 0;
						if(k[i-wI-1] == 1)
							wp++;
						else bp++;
						if(k[i-wI] == 1)
							wp++;
						else bp++;
						if(k[i-1] == 1)
							wp++;
						else bp++;
						if(wp > bp) k[i] = 1;
						if(wp < bp)	k[i] = 0;
					}
				}
				for(int i = 0; i< wI*hI; i++)
					if(k[i] == 0) k[i] = 1;
					else k[i] = 0;
	return k;

}