velocity
========

Extract velocity depth profiles from data file.  

Usage
========

<pre>
usage: velocity.py [-h] [-r LAT LON RADIUS] DATAFILE

Extract velocity/depth profiles from input data file.

positional arguments:
  DATAFILE              Data file from which to extract velocity/depth
                        profiles

optional arguments:
  -h, --help            show this help message and exit
  -r LAT LON RADIUS, --radius LAT LON RADIUS
                        lat,lon and radius (in degrees) for profile search

To select all of the profiles in a 0.25 degree radius around UC Berkeley stadium:

./velocity.py -r 37.869045 -122.264826 0.25 gscnew.txt 
Profile uid: 2206 (37.87N, -122.26E)
Vp	Vs	H	Depth
nan	3.6	17.0	0.0	
nan	4.5	0.0	17.0	


Profile uid: 2248 (37.83N, -122.45E)
Vp	Vs	H	Depth
6.0	nan	20.0	0.0	
8.0	nan	0.0	20.0	


Profile uid: 2257 (37.72N, -122.07E)
Vp	Vs	H	Depth
6.0	nan	28.0	0.0	
8.0	nan	0.0	28.0

</pre>