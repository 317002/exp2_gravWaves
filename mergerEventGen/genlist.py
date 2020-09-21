'''This generates the text file that specifies the merger events for the
processing script.

output file format:
event1ID detectorid1 detectorid2 detectorid3 ...
event2ID detectorid1 detectorid2 detectorid3 ...
event3ID detectorid1 detectorid2 detectorid3 ...
.
.
.

'''
from collections import namedtuple

####Prodject Specific Imports

#used to find events from a catalog
from gwosc import catalog
#used for determing measurment infromation for merger events
from gwosc import datasets


####Constant Perameters####

#This id string points to a collection of observed merger events that the powers
#that be are confident are infact merger events
eventCatalog = 'GWTC-1-confident'



def main():
    #init global constants
    global eventCatalog
    #gets a list of merger events that are part of the eventCatalog var
    events = catalog.events(eventCatalog)
    #making dic object where the event id points to a list of the detectors
    #that measured it
    events = {str(c):datasets.event_detectors(c) for c in events}
    with open('listOfMergerEvents.txt','a+') as doc:
        doc.write(eventCatalog + '\n')
        #for each merger event id
        for event in events:
            detectorString = ''
            #for each detector for each merger event
            for d in events[event]:
                detectorString += ' {}'.format(d)
            #removes the white space in the begining
            detectorString = detectorString[1:]

            doc.write('{} {}\n'.format(event,detectorString))





if __name__ == '__main__':
    main()
