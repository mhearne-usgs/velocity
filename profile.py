#!/usr/bin/env python

#stdlib libraries
import time
import datetime
import argparse
import sys

#third-party libraries
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class velocityProfile:
    '''
    class velocityProfile is ment to handle a single velocity profile from the GSN seismic velocity database.
    and offers depth continuous querying of velocity as well as plotting a stepped V-z plot.
    '''
    def __init__(self, Vp, Vs, h, depth, lat=0, lon=0, uid=False):
        '''
        build-in function __init__ stores input params as numpy.arrays and defines the stepped velocity function

        if there's no velocity information it will be set to numpy.nan

        :param Vp: p-wave velocity as float list in km/s, if n/a .00, type list or numpy.array
        :param Vp: p-wave velocity as float list in km/s, type list or numpy.array
        :param h: layer thickness, type list or numpy.array
        :param lat: profiles latitude, type float
        :param lon: profile longitude, type float
        :param uid: unique id from GSN list, type int
        '''
        self.vp = np.array(Vp)
        self.vs = np.array(Vs)
        self.h = np.array(h)
        self.depth = np.array(depth)
        self.lat = np.float(lat)
        self.lon = np.float(lon)
        self.uid = int(uid)


        self.dstep = np.zeros(self.depth.size*2)
        for i in range(self.depth.size):
            self.dstep[2*i] = self.depth[i]
            if not i == self.depth.size-1:
                self.dstep[2*i+1] = self.depth[i+1]
            else:
                self.dstep[2*i+1] = self.depth[i]+2

        self.vpstep = self.__velocsteps(self.vp)
        self.vsstep = self.__velocsteps(self.vs)

        self.vsstep[self.vsstep==0] = np.nan
        self.vpstep[self.vpstep==0] = np.nan

        self.vs[self.vs==0] = np.nan
        self.vp[self.vp==0] = np.nan
        #self.h[self.h==0] = np.nan

    def __velocsteps(self, v):
        '''
        private function __velocsteps returns a stepped velocity array for :param v:
        used for plotting velocities

        :param v: velocity list, type numpy.array
        :return: stepped velocities, type numpy.array
        '''
        vstep = np.empty(v.size*2)
        for i in range(v.size):
            vstep[2*i] = v[i]
            vstep[2*i+1] = v[i]
        return vstep

    def veloc_at_depth(self, depth, phase="p"):
        '''
        function veloc_at_depth returns a continuous velocity function over depth

        :param depth: list or numpy.array of depths
        :param phase: wave phase, P or S, type string
        :return: array of depth corresponding velocities, type numpy.array
        '''
        if phase == "p":
            vphase = self.vp
        elif phase == "s" and np.any(self.vs):
            vphase = self.vs
        else:
            return np.nan

        ret_v = np.empty(len(depth))
        ret_v[:] = np.nan
        for ri in range(len(depth)):
            for i in range(self.depth.size-1):
                if depth[ri] >= self.depth[i] and depth[ri] <= self.depth[i+1]:
                    ret_v[ri] = vphase[i]
                elif depth[ri] >= self.depth[i+1] and depth[ri] <= (self.depth[i+1] + 2):
                    # Halfspace
                    ret_v[ri] = vphase[i+1]

        return ret_v

    def plot(self):
        '''
        function plot shows the velocity - depth function through matplotlib
        '''
        plt.plot(self.vpstep,-self.dstep, marker='o', color='blue', linestyle= '-')
        if np.any(self.vs):
            plt.plot(self.vsstep,-self.dstep ,marker='o', color='red', linestyle= '-')
        plt.title('Crustal Velocity Structure at ' + str(self.lat) + " / " + str(self.lon))
        plt.xlabel('Velocity [km/s]')
        plt.ylabel('Depth [km]')
        plt.show()

    def __str__(self):
        '''
            built-in function __str__ returns a string to use for print

            :return: multiline string containing all velocity, thickness and geographical information, type string
        '''
        output = "Profile uid: " + str(self.uid) + " ("+ str(self.lat) + "N, " + str(self.lon) + "E)"
        output += "\nVp\tVs\tH\tDepth\n"
        for l in np.arange(len(self.h)):
            output += str(self.vp[l]) + "\t" + str(self.vs[l]) + "\t" + str(self.h[l]) + "\t" + str(self.depth[l]) + "\t\n"
        return output


class containerProfiles:
    '''
    class containerProfiles extents class velocityProfile as it gathers velocity profiles
    and provides functions for spatial selection, querying and processing of the data.
    '''
    def __init__(self):
        '''
        built-in function __init__
        '''
        self.profiles = []
        self.selection = False

        '''
        init list for selection lats and lons, velocities, thicknesses and depths
        '''
        self.lats = False; self.lons = False
        self.vp_array = False; self.vs_array = False; self.h_array = False; self.d_array = False

    def __len__(self):
        '''
        built-in function __len__

        :return: number of profiles, type int
        '''
        return len(self.profiles)

    def __setitem__(self, key, value):
        '''
        built-in function __setitem__
            ? instead of False an exception might has to be triggered

        :param key: array key, type int
        :param value: new value, type velocityProfile()
        '''
        if not value.__class__.__name__ == 'velocityProfile':
            return False
        self.profiles[key] = value

    def __delitem__(self, key):
        '''
        built-in function __ del__

        :param key: array key, type int
        '''
        self.profiles.remove(key)

    def __getitem__(self, key):
        '''
        built-in function __getitem__

        :param key: array key, type int
        '''
        return self.profiles[key]

    def __str__(self):
        '''
        built-in function __str__ to use for print

        :return: multiline string stats about the container, type string
        '''
        output = "Container contains " + str(len(self.profiles)) + " velocity profiles:\n\n"
        output += "uid\tlat\tlon\tmax depth\n"
        #for profile in self.profiles:
        #    output += str(profile.uid) + '\t'+ str(profile.lat) + '\t'+ str(profile.lon) + '\t'+ str(max(profile.depth)) + '\n'
        return output

    def append(self, value):
        '''
        function append handles appending profiles.
            ? instead of False an exception might has to be triggered

        :param value: valocity profile to append to self.profiles, type velocityProfile()
        '''
        if not value.__class__.__name__ == 'velocityProfile':
            return False
        self.profiles.append(value)

    ## Container functions

    def plot_v(self, selection=[]):
        '''
        function plot_v plots velocity depth funtions for :param selection: or the whole dataset

        :param selection: list of selected profiles, type numpy.array or list()
        '''
        if not np.any(selection):
            if len(self.profiles) < 50:
                selection = range(len(self.profiles))
            else:
                selection = range(50)
        plot_profiles = [ self.profiles[i] for i in selection ]
        # Plot selected profiles
        for profile in plot_profiles:
            plt.plot(profile.vpstep,-profile.dstep, marker='o', color='blue', linestyle= '-')
            if np.any(profile.vs):
                plt.plot(profile.vsstep,-profile.dstep ,marker='o', color='red', linestyle= '-')
        plt.title('Crustal Velocity Structure for ' + str(len(plot_profiles)) + " Profiles")
        plt.xlabel('Velocity [km/s]')
        plt.ylabel('Depth [km]')
        plt.show()

    def init_data(self, selection=[]):
        '''
        function init_data initializes continuous arrays

        :param selection: list of selected profiles, type list or numpy.array
        if selection is blank all profiles will be queried

        :self vp_array: array of all P wave velocities, type numpy.array
        :self vs_array: array of all S wave velocities, type numpy.array
        :self h_array: array of all thicknesses, type numpy.array
        :self d_array: array of all depths, type numpy.array
        :self lats: array of all latitudes, type numpy.array
        :self lons: array of all longitudes, type numpy.array
        '''
        if not np.any(selection):
            selection = range(len(self.profiles))
        selected_profiles = [ self.profiles[i] for i in selection ]

        vp_list = []; vs_list = []; h_list = []; d_list = []; lat_list = []; lon_list = []

        for profile in selected_profiles:
            vp_list.append(profile.vp)
            vs_list.append(profile.vs)
            h_list.append(profile.h)
            d_list.append(profile.depth)
            lat_list.append(profile.lat)
            lon_list.append(profile.lon)

        self.vp_array = np.concatenate(vp_list)
        self.vs_array = np.concatenate(vs_list)
        self.h_array = np.concatenate(h_list)
        self.d_array = np.concatenate(d_list)
        self.lats = np.array(lat_list)
        self.lons = np.array(lon_list)

    def plot_hist(self, range=(6, 7), bins=10, phase="p"):
        '''
        function plot_hist plots a velocity histogram  within :param range: and :param bins:

        ! same could be accomplished querying function init_data and analysing self.vp and self.h

        :param range: frequency range list (min, max), type list
        :param bins: number of bins, type int
        '''
        if not np.any(self.vp_array) or np.any(self.vs_array):
            self.init_data()

        v_data = self.vp_array
        plt.hist(v_data, range=range, bins=bins, weights=self.h_array)
        plt.show()

    def depth_at_veloc(self, v_max, dz=.2, d_max=60, d_min=20, phase="p"):
        '''
        function depth_at_veloc returns the last depth :param v_max: has not been exceeded.

        :param v_max: maximal velocity, type float
        :param dz: depth is sampled in dz steps, type float
        :param d_max: maximum depth, type int
        :param d_min: minimum depth, type int
        :param phase: phase to query for, type string NOT YET IMPLEMENTED!!!

        :return: list of lat, lon, depth and uid when v_max is exceeded, type list(np.array)
        '''
        if not self.selection:
            selection = range(len(self.profiles))
        else: selection = self.selection

        lats = np.empty(len(selection))
        lons = np.empty(len(selection))
        depths = np.empty(len(selection))
        uids = np.empty(len(selection), int)
        kicked = 0
        lats[:] = np.nan; lons[:] = np.nan; depths[:] = np.nan
        d_velocities = np.linspace(d_min, d_max, (d_max-d_min)/dz)

        for i in range(len(selection)):
            p = selection[i]
            p_velocities = self.profiles[p].veloc_at_depth(d_velocities)
            for d in range(len(d_velocities)):
                if p_velocities[d] > v_max:
                    lats[i]= self.profiles[p].lat
                    lons[i] = self.profiles[p].lon
                    uids[i] = self.profiles[p].uid
                    if d-1 > 0:
                        depths[i] = d_velocities[d-1]
                    else:
                        depths[i] = d_velocities[d]
                    break
        print "kicked", kicked
        return lats[depths>=d_min], lons[depths>=d_min], depths[depths>=d_min], uids[depths>=d_min]


    ## Selection Functions
    def select_region(self, west, east, south, north):
        '''
        function select_region selects a region by geographic coordinates :param west, east, south, north:

        :param west: west edge of region, type float
        :param east: east edge of region, type float
        :param south: south edge of region, type float
        :param north: north edge of region, type float

        :return, self selection: array containing all profiles keys within desired region, type numpy.array
        '''
        # Select Region by lat and lon
        #
        self.vp_array = False; self.vs_array = False; self.h_array = False; self.d_array = False
        self.selection = []

        for i in range(len(self.profiles)):
            if self.profiles[i].lon >= west and self.profiles[i].lon <= east \
            and self.profiles[i].lat <= north and self.profiles[i].lat >= south:
                self.selection.append(i)
        return np.array(self.selection)

    def select_plotv(self):
        '''
        function plot_v plots velocity depth funtions for :self selection:
        '''
        self.plot_v(selection=self.selection)

    def select_data(self):
        '''
        function select_data initializes the continuous arrays

        see function init_data
        '''
        self.init_data(selection=self.selection)

    def select_poly(self, poly):
        '''
        function select_poly determines if a profile is inside a given polygon or not
        polygon is a list of (x,y) pairs

        The algorithm is called the "Ray Casting Method"

        :param poly: list of x, y pairs, type list(numpy.array(2))

        :return, self selection: array containing all profiles keys within desired region, type numpy.array
        '''
        selection = []

        n = len(poly)

        for n in range(len(self.profiles)):
            x = self.profiles[n].lat
            y = self.profiles[n].lon

            inside = False
            p1x,p1y = poly[0]
            for i in range(n+1):
                p2x,p2y = poly[i % n]
                if y > min(p1y,p2y):
                    if y <= max(p1y,p2y):
                        if x <= max(p1x,p2x):
                            if p1y != p2y:
                                xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                            if p1x == p2x or x <= xints:
                                inside = not inside
                p1x,p1y = p2x,p2y
            if inside:
                selection.append(n)

        self.selection = np.array(selection)

    def select_depth(self, depth=30):
        '''
        function select_depth selects only that profiles deeper than :param depth:
        either for the whole dataset or, if set, for the selected datasets

        :param depth: minimum depth for the profiles, type float

        :return, self selection: array containing all profiles keys deeper than desired depth, type numpy.array
        '''
        if not self.selection:
            selection = range(len(self.profiles))
        else: selection = self.selection

        selected = []
        for profile in selection:
            if max(self.profiles[profile].depth) >= depth:
                selected.append(profile)

        self.selection = selected
        return self.selection

    def select_location(self, lat, lon, radius=10):
        '''
        function select_location selects profiles at :param lat, lon: within a :param radius:

        :param lat: latitude, type float
        :param lon: longitude, type float
        :param radius: radius surrounding lat, lon, type float

        :return, self selection: array containing all profiles within radius of location lat, lon, type numpy.array
        '''
        selected = []
        for p in range(len(self.profiles)):
            if np.sqrt((lat - self.profiles[p].lat)**2 + (lon - self.profiles[p].lon)**2) <= radius:
                selected.append(p)

        self.selection = selected
        return self.selection


    def select_delete(self):
        '''
        function select_delete deletes the selections made

        :return: True, type bool
        '''
        self.selection = False
        return True


def readProfiles(database_file):
    '''
    function readProfiles reads in the the GSN databasefile and puts it in containerProfiles

    !!TO BE INCLUDED INTO containerProfiles!!


    File format:

       uid  lat/lon  vp    vs    hc     depth
        2   29.76N   2.30   .00   2.00    .00  s  25.70   .10    .00  NAC-CO   5 U
            96.31W   3.94   .00   5.30   2.00  s  33.00   MCz  39.00  61C.3    EXC
                     5.38   .00  12.50   7.30  c
                     6.92   .00  13.20  19.80  c
                     8.18   .00    .00  33.00  m

        3   34.35N   3.00   .00   3.00    .00  s  35.00  1.60    .00  NAC-BR   4 R
           117.83W   6.30   .00  16.50   3.00     38.00   MCz  55.00  63R.1    ORO
                     7.00   .00  18.50  19.50
                     7.80   .00    .00  38.00  m


    :param database_file: path to database file, type string

    :return: class containerProfiles() containing all entries from database, type containerProfiles()
    '''

    with open(database_file, 'r') as database:
        vp=[]; vs=[]; h=[]; depth=[]; readline=1; profiles = containerProfiles()
        db = database.readlines()
        for dbline in db:
            if dbline.isspace():
                # create velocity profile from read lines
                if not len(depth) == 0:
                    profiles.append(velocityProfile(vp, vs, h, depth, lat, lon, uid))
                if not len(vp) == len(h):
                    print uid, lat, lon, vp, vs, h, depth
                vp=[]; vs=[]; h=[]; depth=[]; readline=1
            else:
                # read in profile data
                if readline == 1:
                    uid = int(dbline[0:6])
                    lat = float(dbline[8:13])
                    if dbline[13] == "S":
                        lat = -lat
                if readline == 2:
                    lon = float(dbline[7:13])
                    if dbline[13] == "W":
                        lon = -lon

                try:
                    vp.append(float(dbline[17:21]))
                    vs.append(float(dbline[23:27]))
                    h.append(float(dbline[28:34]))
                    depth.append(float(dbline[35:41]))
                except:
                    pass

                readline += 1
        # Append last profile
        profiles.append(velocityProfile(vp, vs, h, depth, lat, lon, uid))
    return profiles

def main(args,parser):
    if not args.datafile:
        print 'Missing required data file name.  Exiting.'
        parser.print_help()
        sys.exit(0)
    if len(args.radius) < 3:
        print 'To do a radius search, you must specify a lat,lon and radius.'
        parser.print_help()
        sys.exit(0)
    worldProfiles = readProfiles(args.datafile)
    worldProfiles.select_location(args.radius[0],args.radius[1],args.radius[2])
    worldProfiles.plot_v(selection=worldProfiles.selection)
    for p in worldProfiles.selection:
        print worldProfiles.profiles[p]
        print
    if args.doPlot:
        fname = 'velocity'+datetime.datetime.now().strftime('%Y%m%d%H%M%S')+'.png'
        plt.savefig(fname)
        sys.stderr.write('Saved plot to %s\n' % fname)

if __name__ == '__main__':
    usage = 'Extract velocity/depth profiles from input data file.'
    cmdparser = argparse.ArgumentParser(description=usage)
    cmdparser.add_argument('datafile', metavar='DATAFILE',
                           help='Data file from which to extract velocity/depth profiles')
    cmdparser.add_argument("-r", "--radius", dest="radius",nargs=3,
                           help="lat,lon and radius (in degrees) for profile search", metavar=('LAT','LON','RADIUS'),type=float)
    cmdparser.add_argument("-p", "--plot", dest="doPlot", default=False,action='store_true',
                           help="Create a PNG plot of results (filename velocity_YYYYMMDDHHMMSS.png')")
    cmdargs = cmdparser.parse_args()
    main(cmdargs,cmdparser)    


# '''
# Database is read in and crustal depth through velocities are queried
# '''
# t = time.time()
# worldProfiles = readProfiles("gscnew.txt")
# print "loaded..."
# print str(time.time() - t)
# # selecting US american region
# worldProfiles.select_region(-130, -60, 20, 55)
# t = time.time()
# lats, lons, depths, uids = worldProfiles.depth_at_veloc(6.35, d_min=0, dz=.1, d_max=30)
# np.savetxt("Vp_6.35.txt", np.array([lons, lats, depths, uids]).T)
# print "Vp_6.35.txt in ", str(time.time() - t)
# t = time.time()
# lats, lons, depths, uids = worldProfiles.depth_at_veloc(6.75, d_min=0, d_max=65)
# np.savetxt("Vp_6.75.txt", np.array([lons, lats, depths, uids]).T)
# print "Vp_6.75.txt in ", str(time.time() - t)
# t = time.time()
# lats, lons, depths, uids = worldProfiles.depth_at_veloc(7.6, d_min=0, d_max=65)
# np.savetxt("Vp_7.6.txt", np.array([lons, lats, depths, uids]).T)
# print "Vp_7.6.txt in ", str(time.time() - t)

