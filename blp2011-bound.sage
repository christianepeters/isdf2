# author: Christiane Peters, http://christianepeters.wordpress.com
# April 2013
# run using 'sage blp2011-bound.sage'

def log2(x):
  return log(x)/log(2)

def myfactorial(x):
  return gamma(x+1)

def mybinomial(x,y):
  return myfactorial(x)/(myfactorial(y)*myfactorial(x-y))

def bound(n,k,w):
  min=[10000,-1]
  for p in range(0,w):
    wf=log2(mybinomial(n,w))-1\
	   -log2(mybinomial(n-k,w-p))\
	   -log2(sqrt(mybinomial(k,p)));
    if wf < min[0]:
      min[0]=RDF(wf); min[1]=p;
  return min;


n=6624;k=5129;w=117
print n,k,w,bound(n,k,w)
n=1024;k=524;w=50;
print n,k,w,bound(n,k,w)
