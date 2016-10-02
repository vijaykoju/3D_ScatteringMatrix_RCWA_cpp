/* va_arg example */
#include <iostream>      /* printf */
#include <cstdarg>     /* va_list, va_start, va_arg, va_end */
using namespace std;

typedef struct input{
	int n;
	int nn;
} in;

int FindMaxStruct(in s_in, ...)
{
	in val, largest;
	int max;
	va_list vl;
	va_start(vl, s_in);
	val = va_arg(vl, in);
	max = val.nn;
	for (int i=1; i<s_in.n; i++)
	{
		val = va_arg(vl, in);
		max = (max>val.nn)?max:val.nn;	
	}
	va_end(vl);
	return max;
}


int FindMax (int n, ...)
{
  int i,val,largest;
  va_list vl;
  va_start(vl,n);
  largest=va_arg(vl,int);
  for (i=1;i<n;i++)
  {
    val=va_arg(vl,int);
    largest=(largest>val)?largest:val;
  }
  va_end(vl);
  return largest;
}

int main ()
{
  int m,n;
  m= FindMax (7,702,422,631,834,892,104,772);
	n = FindMaxStruct((in){7,702}, (in){2,422}, (in){2,631}, (in){2,834}, (in){2,892}, (in){2,104}, (in){2,772});
	cout << "The largest value is: " << m << endl;
	cout << "The largest value is: " << n << endl;
  //printf ("The largest value is: %d\n",m);
  return 0;
}
