#!/usr/bin/env python

'''An experimental module containing the Make_Journey class
which plots a world map and a great circle journeys on
that map. J Line 2018'''

from __future__ import division,print_function
from numpy import *
import matplotlib
matplotlib.use('agg')
matplotlib.rcParams.update({'font.size': 18})
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from subprocess import call
import cv2
from datetime import datetime

##Here are a number of latitide/longitude pairs

newyork = [40.78,-73.98]
melbourne = [-37.813611, 144.963056]
london = [51.53,0.08]
seattle = [47.609722, -122.333056]
losangeles = [34.05, -118.25]
riodj = [-22.908333, -43.196389]
singapore = [1.283333, 103.833333]
auckland = [-36.840556, 174.74]
sydney = [-33.865, 151.209444]
bangkok = [13.75, 100.466667]
cairns = [-16.092222, 145.775278]
darwin = [-12.45, 130.833333]
dubai = [24.95, 55.333333]
baries = [-34.603333, -58.381667]
igazu = [-25.686667, -54.444722]
salta = [-24.783333, -65.416667]
lapaz = [-16.5, -68.15]
lima = [-12.043333, -77.028333]
toronto = [43.7, -79.4]
tenerife = [28.268611, -16.605556]
heraklion = [35.339722, 25.180278]
split = [43.538889, 16.298056]
deuxalpes = [45.007222, 6.121667]
valdisere = [45.4506, 6.9781]
bratislava = [48.143889, 17.109722]
malta = [35.883333, 14.5]
tallinn = [59.437222, 24.745278]
uluru = [-25.345, 131.036111]
brisbane = [-27.466667, 153.033333]
longreach = [-23.45, 144.25]
adelaide = [-34.929, 138.601]
dunedin = [-45.866667, 170.5]
sanfran = [37.783333, -122.416667]
vancouver = [49.25, -123.1]
perth = [-31.952222, 115.858889]
MWA = [-26.703319, 116.670815]
detroit = [42.2125, -83.353333]
pennstate = [40.849167, -77.848611]
changi = [1.359211, 103.989306]
calgary = [51.05, -114.066667]
doha = [25.286667, 51.533333]
pheonix = [33.434167, -112.011667]
medhat = [50.041667, -110.6775]
vienna = [48.2, 16.366667]
innsbruck = [47.260278, 11.343889]
munich = [48.133333, 11.566667]
thessaloniki = [40.519722, 22.970833]
athens = [37.936389, 23.947222]
mumbai = [18.975, 72.825833]
GMRT = [19.096517, 74.049742]
houston = [29.984444, -95.341389]
alburquerque = [35.039333, -106.610778]
denver = [39.861667, -104.673056]
chicago = [41.978611, -87.904722]
guangzhou = [23.3925, 113.298889]
klumpar = [2.743333, 101.698056]
taiwan = [25.033333, 121.633333]
hongkong = [22.3, 114.2]
amsterdam = [52.3702, 4.8952]
venice = [45.4408, 12.3155]
tokyo = [35.6895, 139.6917]
kumamoto = [32.8031, 130.7079]
shikoku = [33.7432, 133.6375]
muscat = [23.5859, 58.4059]
kathmandu = [27.7172, 85.3240]
klumpuar = [3.1390, 101.6869]
helsinki = [60.1699, 24.9384]
jakarta = [6.1751, 106.8650]


class Make_Journey(object):
    '''This Class uses Basemap from mpl_toolkits to plot a world map.
    The methods within allow a user to plot a custom journey on that map
    The ouputs are save in a ./plot folder, which is created if necessary.
    A plot is created for each step, which can be used to create a movie'''
        
    def _save(self,leg_index):
        '''Saves the current instance of our figure - clears
        previous text and adds updated text like current location
        and distance travelled'''
        
        ##Clear previous text about the figure
        for txt in self.fig.texts:
            txt.set_visible(False)

        ##Add the current distance travelled
        self.fig.text(0.09, 0.1, 'Distance:\n{:,d}$\,$km'.format(int(round(self._distance))),
        verticalalignment='top', horizontalalignment='center',fontsize=28)
        
        ##If information available, add the current location label
        ##top left, and the date bottom right
        if self.from_file:
            name = self.names[leg_index]
            
            try:
                label = self.label_dict[name]
                self.fig.text(0.82, 0.1,'Location:\n%s' %label,verticalalignment='top', horizontalalignment='left',fontsize=28)
            except:
                pass
            
            if self.dates:
                month,year = self.dates[leg_index].strftime('%B %Y').split()
                self.fig.text(0.05, 0.95,'%s\n%s' %(month,year),verticalalignment='top', horizontalalignment='left',fontsize=28)

        ##Draw the current canvas - usually gets done with either plt.show()
        ##or fig.savefig, but we're saving straight into the video object
        self.fig.canvas.draw_idle()
        
        ##Convert plot into a numpy array
        data = fromstring(self.fig.canvas.tostring_rgb(), dtype=uint8)
        ##Reshape that array into 2D x red, green, blue colours (hence adding a 3rd dimension of 3)
        frame = data.reshape(self.fig.canvas.get_width_height()[::-1] + (3,))
        
        ##3rd dimension comes out backwards, makes colour weird - reverse it
        self.video.write(frame[:,:,::-1])
        
        ##If the user wants to save every frame, save them in the 'plots' dir
        if self.save_all_figs:
            self.fig.savefig('./plots/flights_%06d.png' %(self.plot_count))
            
    def _read_lonlat(self,filename):
        '''Reads a text file (named filename) which should contain lon/lat coords
        Turns them into dictionaries to be used in the class'''
        
        latlon_dict = {}
        label_dict = {}
        
        lines = [line for line in open(filename,'r').read().split('\n') if line!='' or '#' in line]
        
        for line in lines:
            split = line.split()
            
            if len(split) == 3:
                name,lat,lon = split
                latlon_dict[name] = [float(lat),float(lon)]
            else:
                name,lat,lon = split[:3]
                label = split[3]
                for extra in split[4:]: label += ' %s' %extra
                
                latlon_dict[name] = [float(lat),float(lon)]
                label_dict[name] = label
                
        self.latlon_dict = latlon_dict
        self.label_dict = label_dict
        
    def _create_journey(self,filename):
        '''Using the dictionaries created by self._read_lonlat,
        auto creates the legs necessary to create the journey
        (Should) automatically detect whether dates are available
        and add that information accordingly'''
        legs = loadtxt(filename,dtype=str)
        
        if len(legs.shape) == 1:
            self.names = [name for name in legs]
            self.multi_list = [self.latlon_dict[name] for name in legs]
        elif len(legs.shape) == 2:
            self.names = [name for name in legs[:,0]]
            self.multi_list = [self.latlon_dict[name] for name in legs[:,0]]
            self.dates = [datetime.strptime(date.replace('/',' '),'%m %Y') for date in legs[:,1]]
            
    def __init__(self,trav_c='r',done_c='#FFA500',marker_c='y',frame_rate=1,fig_height=14,fig_width=24,save_name='flights',save_all_figs=False,from_file=False,locations=False):
        '''
        Sets up the map with default colours and sizes. If from_file and
        locations are set, creates a journey movie using those
        
        trav_c (default='r') = colour of great circle during travel
        done_c (default='#FFA500') = colour of great circle after travel
        marker_c (default='y') = colour of location marker
        frame_rate (default='1') = number of frames per second for vidoes
        fig_height (default=14) = height of figure (inches)
        fig_width (default=24) = width of figure (inches)
        save_name (default='flights') = prefix name for the outputs 
        save_all_figs (default=False) = If True, save every frame update to a new figure
        from_file (default=False) = 
        
        '''
        ##Setup initial conditions
        self.plot_count = 0
        self.trav_c = trav_c
        self.done_c = done_c
        self.marker_c = marker_c
        ##Set the number of frames per second for the video
        self.frame_rate = frame_rate
        self.save_name = save_name
        ##Initialise the distance counter at zero
        self._distance = 0
        ##If True, each frame is saved to a Figure
        self.save_all_figs = save_all_figs
        ##If True, make journey from a file
        self.from_file = from_file
        self.dates = False
        
        ##Create a figure and axis for plot
        self.fig = plt.figure(figsize=(fig_width,fig_height))
        self.ax = self.fig.add_axes([0.04,0.04,0.92,0.92])
        
        ##Create the basemap figure
        self.m = Basemap(lat_0=0.0,lon_0=0.0,projection='kav7', resolution ='c',ax=self.ax)
        self.m.drawcoastlines()
        self.m.drawcountries()
        self.m.drawmapboundary(fill_color='c')
        self.m.fillcontinents(color='#B3A580',lake_color='c')
        self.m.drawparallels(array([-50.0,-25.0,0.0,25.0,50.0,75.0]),labels=[1,1,0,0])
        self.m.drawmeridians(arange(-180,210,45),labels=[0,0,1,1])
        
        fourcc = cv2.VideoWriter_fourcc(*'mp4v')
        
        self.video = cv2.VideoWriter('%s_flight_movie.mp4' %self.save_name, fourcc, self.frame_rate, (self.fig.canvas.get_width_height()))
        
        if self.from_file:
            ##Makes a dictionary containing
            self._read_lonlat(locations)
            self._create_journey(from_file)
            
        ##Save the blank map plot
        self._save(leg_index=0)
            
        if self.from_file:
            self.multi_leg(self.multi_list)
            self.finish_journey()
        
        
    def _haversine(self,lon_1=None,lat_1=None,lon_2=None,lat_2=None):
        '''Calculates the distance in kms between two lat/lon
        points (as the crow flies) using the Haversine formula
        https://en.wikipedia.org/wiki/Haversine_formula
        args: lon_1,lat_1,lon_2,lat_2 in degrees
        '''
        ##Convert degree to radian
        D2R = pi / 180.
        ##Radius of the Earth in km
        R = 6371.

        ##Set up coords in radian
        lamb1 = lon_1 * D2R
        lamb2 = lon_2 * D2R
        phi1 = lat_1 * D2R
        phi2 = lat_2 * D2R
        phi_diff = phi2 - phi1
        lamb_diff = lamb2 - lamb1
        
        return 2*R*arcsin(sqrt(sin(phi_diff/2.)*sin(phi_diff/2.) + cos(phi1)*cos(phi2)*sin(lamb_diff/2.)*sin(lamb_diff/2.)))
        
    def _draw_line(self,coords_1=None,coords_2=None,colour=None,plot_point=False,save_plot=False,calc_distance=True,leg_index=None):
        '''This draws a great circle from coords_1 to coords_2
        coords_1 = list(longitude1, latitude1) in deg
        coords_2 = list(longitude2, latitude2) in deg
        '''
        lat_1,lon_1 = coords_1
        lat_2,lon_2 = coords_2
        line, = self.m.drawgreatcircle(lon_1,lat_1,lon_2,lat_2,linewidth=2.0,color=colour,del_s=100)
        #Line below might be useful in the future
        #line_coords = line.get_path().vertices
        
        if calc_distance:
            ##calcualte the distance between the two points
            self._distance += self._haversine(lon_1,lat_1,lon_2,lat_2)
        
        if plot_point:
            self.m.scatter(lon_2,lat_2,50,marker='o',edgecolor='k',color=self.marker_c,linewidth=0.5,latlon=True,zorder=100)
            
        self.last_coords = [coords_1,coords_2]
        #self.last_leg_index
        
        if save_plot:
            self.plot_count += 1
            print('Doing leg',self.plot_count)
            self._save(leg_index)
            
    def return_journey(self,coords_1=None,coords_2=None):
        '''Adds a return journey between ther lists coords_1 and coords_2
        Args:
        coords_1 = list(longitude1, latitude1) in deg
        coords_2 = list(longitude2, latitude2) in deg
        '''
        
        self._draw_line(coords_1=coords_1,coords_2=coords_2,colour=self.trav_c,plot_point=True,save_plot=True)
        ##Although we are plotting a 'done' leg here, we are travelling back from the desination,
        ##so we don't set calc_distance=False because we want to include this disance
        self._draw_line(coords_1=coords_2,coords_2=coords_1,colour=self.done_c,save_plot=True)
            
    def multi_leg(self,coords=None):
        self._draw_line(coords_1=coords[0],coords_2=coords[1],colour=self.trav_c,plot_point=True,save_plot=True,leg_index=1)
        for i in arange(len(coords)-2):
            self._draw_line(coords_1=coords[i],coords_2=coords[i+1],colour=self.done_c,calc_distance=False,leg_index=i+1)
            
            ##If doing a return journey, end up with two red lines in a row, looks rubbish. Check using
            ##the last set of coords - if doing return, use the colour done_c instead
            if coords[i+2] == self.last_coords[0]:
                self._draw_line(coords_1=coords[i+1],coords_2=coords[i+2],colour=self.done_c,plot_point=True,save_plot=True,leg_index=i+2)
            else:
                self._draw_line(coords_1=coords[i+1],coords_2=coords[i+2],colour=self.trav_c,plot_point=True,save_plot=True,leg_index=i+2)
            
        self._draw_line(coords_1=coords[-2],coords_2=coords[-1],colour=self.done_c,calc_distance=False,leg_index=-1)
        
    def finish_journey(self):
        self._draw_line(coords_1=self.last_coords[0],coords_2=self.last_coords[1],colour=self.done_c,save_plot=True,calc_distance=False,leg_index=-1)
        
        self.fig.savefig('%s_full_journey.png' %self.save_name,bbox_inches='tight')
        

    #
if __name__=='__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Plot')
    parser.add_argument('integers', metavar='N', type=int, nargs='+',
                        help='an integer for the accumulator')
    parser.add_argument('--sum', dest='accumulate', action='store_const',
                        const=sum, default=max,
                        help='sum the integers (default: find the max)')

    args = parser.parse_args()
    
    
    journey = Make_Journey(save_all_figs=False,save_name='example_fromtxt',from_file='example_journey.txt',locations='lonlat_locations.txt')

    