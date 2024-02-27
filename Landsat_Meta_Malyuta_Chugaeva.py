import os
from math import radians, cos, sin, asin, sqrt
from geopy.geocoders import Nominatim

path = (os.path.dirname(os.path.abspath(__file__)))
os.chdir(path)

# Make sure to place the Landsat 8/9 metadata .xml file in the same folder as the script

with open('LC08_L2SP_225084_20240214_20240223_02_T1_MTL.xml') as f: # Enter your metadata file name here
    meta = f.read()
lattitudes = ['<CORNER_UL_LAT_PRODUCT>', '</CORNER_UL_LAT_PRODUCT>', 
              '<CORNER_UR_LAT_PRODUCT>', '</CORNER_UR_LAT_PRODUCT>', 
              '<CORNER_LL_LAT_PRODUCT>', '</CORNER_LL_LAT_PRODUCT>', 
              '<CORNER_LR_LAT_PRODUCT>', '</CORNER_LR_LAT_PRODUCT>']   # list of indicators, between which the latitudes are located 
lattitudesvalues = []
for coordinate in lattitudes:
    if '/' not in coordinate:
        index1 = meta.find(coordinate) + len(coordinate)   # index1 -- head index of the required value
        continue
    index2 = meta.find(coordinate)        # index2 -- final index of the required value
    lattitudesvalues.append(float(meta[index1:index2]))      # list of latitudes

longtitudes = ['<CORNER_UL_LON_PRODUCT>', '</CORNER_UL_LON_PRODUCT>', 
               '<CORNER_UR_LON_PRODUCT>', '</CORNER_UR_LON_PRODUCT>', 
               '<CORNER_LL_LON_PRODUCT>', '</CORNER_LL_LON_PRODUCT>', 
               '<CORNER_LR_LON_PRODUCT>', '</CORNER_LR_LON_PRODUCT>']    # list of indicators, between which the longtitudes are located 
longtitudesvalues = []
for coordinate in longtitudes:
    if '/' not in coordinate:
        index1 = meta.find(coordinate) + len(coordinate)
        continue
    index2 = meta.find(coordinate) 
    longtitudesvalues.append(float(meta[index1:index2]))    # list of longtitudes

def findchars(meta):
    """Return copmosite image qualities such as:
    date, time, sun position, cloud coverage and image quality for different bands: OLI and TIRS

    Keyword arguments:
    meta -- file with metadata
    """
    qualitives = ['<DATE_ACQUIRED>', '</DATE_ACQUIRED>', 
                  '<SCENE_CENTER_TIME>', '</SCENE_CENTER_TIME>', 
                  '<SUN_ELEVATION>', '</SUN_ELEVATION>', 
                  '<CLOUD_COVER>', '</CLOUD_COVER>', 
                  '<IMAGE_QUALITY_OLI>', '</IMAGE_QUALITY_OLI>', 
                  '<IMAGE_QUALITY_TIRS>', '</IMAGE_QUALITY_TIRS>', 
                  '<SPACECRAFT_ID>', '</SPACECRAFT_ID>']   # list of indicators, between which the required values are located 
    qualitivevalues = []
    for quality in qualitives:
        if '/' not in quality:
            indexn1 = meta.find(quality) + len(quality)  # indexn1 -- head index of the required value
            continue
        indexn2 = meta.find(quality)    # indexn2 -- final index of the required value
        qualitivevalues.append(meta[indexn1:indexn2])
    return qualitivevalues  # list of image quality values
quals = findchars(meta) # quals -- final values of each image quality

def haversine(lattitudesvalues, longtitudesvalues, r=6371):

    hvrsines = []
    distances = []
    """Calculates the spatial coverage of the image using their coordinates and haversine formula
    
    Keyword arguments:  
    lattitudesvalues -- list of the lattitudes of 4 points 
    longtitudesvalues -- list of the longtitudes of 4 points 
    r -- Earth radius

    points: [0]*-------*[1]
               | \      |
               |   \    |
               |     \  |
            [2]*-------*[3]
    """
    lattitudesvalues = list(map(radians, lattitudesvalues))    # converts coordinates from degrees to radians
    longtitudesvalues = list(map(radians, longtitudesvalues))

    dlat1 = lattitudesvalues[0] - lattitudesvalues[2]    # calculates the distance between the lattitudes of two points
    dlong1 = longtitudesvalues[0] - longtitudesvalues[2]     # calculates the distance between the longtitudes of two points
    dlat2 = lattitudesvalues[1] - lattitudesvalues[0] 
    dlong2 = longtitudesvalues[1] - longtitudesvalues[0]
    dlat3 = lattitudesvalues[3] - lattitudesvalues[2] 
    dlong3 = longtitudesvalues[3] - longtitudesvalues[2]
    dlat4 = lattitudesvalues[1] - lattitudesvalues[3] 
    dlong4 = longtitudesvalues[1] - longtitudesvalues[3]
    dlat5 = lattitudesvalues[0] - lattitudesvalues[3] 
    dlong5 = longtitudesvalues[0] - longtitudesvalues[3]

    hvrsines.append(sin(dlat1/2)**2 + cos(lattitudesvalues[0]) * cos(lattitudesvalues[2]) * sin(dlong1/2)**2) # substitutes values into the formula to calculate haversine
    hvrsines.append(sin(dlat3/2)**2 + cos(lattitudesvalues[3]) * cos(lattitudesvalues[2]) * sin(dlong3/2)**2)
    hvrsines.append(sin(dlat2/2)**2 + cos(lattitudesvalues[1]) * cos(lattitudesvalues[0]) * sin(dlong2/2)**2)
    hvrsines.append(sin(dlat4/2)**2 + cos(lattitudesvalues[1]) * cos(lattitudesvalues[3]) * sin(dlong4/2)**2)
    hvrsines.append(sin(dlat5/2)**2 + cos(lattitudesvalues[0]) * cos(lattitudesvalues[3]) * sin(dlong5/2)**2) # diagonal

    for value in hvrsines:
        dist = r * 2 * asin(sqrt(value))# calculates the distance between two points; the last iteration calculates the diagonal
        distances.append(dist)
    
    poluperimeter1 = (distances[0] + distances[1] + distances[4]) / 2 # calculates semiperimeter of the first triangle
    poluperimeter2 = (distances[2] + distances[3] + distances[4]) / 2

    square1 = sqrt(poluperimeter1 * (poluperimeter1 - distances[0]) * (poluperimeter1 - distances[1]) * (poluperimeter1 - distances[4])) # calculates the square of the first triangle using Heron's formula
    square2 = sqrt(poluperimeter2 * (poluperimeter2 - distances[2]) * (poluperimeter2 - distances[3]) * (poluperimeter2 - distances[4]))

    return (square1 + square2) # returns spatial coverage of the image


def center(lattitudesvalues, longtitudesvalues):
    """Calculates the center of image

    Keyword arguments:
    lattitudesvalues -- list of the lattitudes of 4 points 
    longtitudesvalues -- list of the longtitudes of 4 points 
    """
    meancenter = []
    meancenter.extend((sum(lattitudesvalues)/4, sum(longtitudesvalues)/4))  # adds average values of latitude and longtitude to the list which is the coordinates of the center
    return meancenter


def hemispherelat(lattitudesvalues):
    """Calculates which hemisphere, northern or southern, the image covers

    Keyword arguments:
    lattitudesvalues -- list of the lattitudes of 4 points 
    """
    for k in range(len(lattitudesvalues)): # Checking whether the image is in the northern or southern hemisphere
        if lattitudesvalues[0] > 0:
            if lattitudesvalues[k] > 0:
                hemisphere_lat = 1
            else:
                hemisphere_lat = 0
                break
        elif lattitudesvalues[0] < 0:
            if lattitudesvalues[k] < 0:
                hemisphere_lat = -1
            else:
                hemisphere_lat = 0
                break
    if hemisphere_lat == 1:
        return 'Nothern'
    elif hemisphere_lat == 0:
        return 'Image covers both Nothern and Southern Hemispheres'
    else:
        return 'Southern'

def hemispherelon(longtitudesvalues):
    """Calculates which hemisphere, eastern or western, the image covers

    Keyword arguments:
    longtitudesvalues -- list of the longtitudes of 4 points 
    """
    for k in range(len(longtitudesvalues)): # Checking whether the image is in the western or eastern hemisphere
        if longtitudesvalues[0] > 0:
            if longtitudesvalues[k] > 0:
                hemisphere_lon = 1
            else:
                hemisphere_lon = 0
                break
        elif longtitudesvalues[0] < 0:
            if longtitudesvalues[k] < 0:
                hemisphere_lon = -1
            else:
                hemisphere_lon = 0
                break
    if hemisphere_lon == 1:
        return 'Eastern'
    elif hemisphere_lon == 0:
        return 'Image covers both Western and Eastern Hemispheres'
    else:
        return 'Western'
        
def geolocate(lattitudesvalues, longtitudesvalues):
    """Identifies in which area the center of the image is located

    Keyword arguments:
    lattitudesvalues -- list of the lattitudes of 4 points 
    longtitudesvalues -- list of the longtitudes of 4 points 
    """
    # initialize Nominatim API
    geolocator = Nominatim(user_agent="Me")
    # Latitude & Longitude input
    Latitude = str(center(lattitudesvalues, longtitudesvalues)[0])
    Longitude = str(center(lattitudesvalues, longtitudesvalues)[1])
    location = geolocator.reverse(Latitude+","+Longitude)
    address = location.raw['address']
    # traverse the data
    state = address.get('state', '')
    country = address.get('country', '')
    return [country, state]
geolocation = geolocate(lattitudesvalues, longtitudesvalues)

def resultexport(lattitudesvalues, longtitudesvalues):
    """Saves the results to a separate new file

    Keyword arguments:
    lattitudesvalues -- list of the lattitudes of 4 points 
    longtitudesvalues -- list of the longtitudes of 4 points 
    """
    with open('result.txt', mode='w+') as res: # name your export file
        res.write("Satellite: " + (str(quals[-1])).capitalize())
        res.write("\n" + 'Image coordinates: ')
        for l in range(len(lattitudesvalues)):    
            res.write('[' + str(lattitudesvalues[l]) + ', ' + str(longtitudesvalues[l]) + ']' + '; ')
        res.write("\n" + 'Center Coordinates: [' + str(center(longtitudesvalues, lattitudesvalues)[1]) + ', ' + str((center(longtitudesvalues, lattitudesvalues)[0])) + ']')
        res.write("\n" + 'N/S Hemisphere: ' + str(hemispherelat(lattitudesvalues)))
        res.write("\n" + 'W/E Hemisphere: ' + str(hemispherelon(longtitudesvalues)))
        res.write("\n" + 'Country Corresponding to the Center: ' + geolocation[0] + ', ' + geolocation[1])
        res.write("\n" + 'Metadata Spatial Coverage (Including nodata values): ' + str(haversine(lattitudesvalues, longtitudesvalues)) + ' km2')
        res.write("\n" + 'Estimated Real Spatial Coverage: ' + str(haversine(lattitudesvalues, longtitudesvalues)*2/3) + ' km2')
        res.write("\n" + 'Image acquisition time: ' + quals[0] + ', ' + quals[1] + ' GMT')
        res.write("\n" + 'Brief qualititive image characteristic:')
        res.write("\n" + '1. Sun Elevation Degree: ' + str(quals[2]))
        if quals[-1] == "LANDSAT_9":
            res.write("\n" + '2. Overall Cloud Cover: ' + str(float(quals[3]) * 100) + '%')
        elif quals[-1] == "LANDSAT_8":
            res.write("\n" + '2. Overall Cloud Cover: ' + str(quals[3]) + '%')
        res.write("\n" + '3. Image Quality (OLI Sensor): ' + str(quals[4]) + '/9')
        res.write("\n" + '4. Image Quality (TIRS Sensor): ' + str(quals[4]) + '/9')
        print('Success! Check the output file in your directory')
resultexport(lattitudesvalues, longtitudesvalues)

