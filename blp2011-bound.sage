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


bound(6624,5129,117)
bound(1024,524,50)
