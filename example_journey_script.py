from numpy import *
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
fig.savefig('',bbox_inches='tight')







journey = Make_Journey(save_all_figs=False,save_name='example_fromtxt',from_file='example_journey.txt',locations='lonlat_locations.txt')





journey = Make_Journey(save_all_figs=False,save_name='example_fromscript')
    #######2008===========================================================================

    ##04 06 2008 - ##11 06 2008
    journey.return_journey(london,split)

    ##02-Sep-2008 - #09-Sep-2008
    journey.return_journey(london,tenerife)

    ###2009===========================================================================
    #09-01-2009 - #16-01-2009
    journey.return_journey(london,deuxalpes)

    ##JUNE 2009
    journey.return_journey(london,bratislava)

    ##04 AUG 2009 - #12 04 AUG 2009
    journey.return_journey(london,heraklion)
        
    ###2010===========================================================================
    ###===============================================================================

    ###JAN 2010
    journey.return_journey(london,valdisere)

    ###JUNE 2010
    journey.return_journey(london,malta)
        
    ##Sun 05 Sep 2010 - Sat, 2 Oct 2010

    journey.multi_leg([london,riodj,igazu,baries,salta,lapaz,lima,toronto,london])
