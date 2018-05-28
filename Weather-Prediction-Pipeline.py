__author__ = 'jammy'

import numpy as np
import matplotlib.pyplot as plt
# Global variables specifying year range
year_min = 1900
year_max = 2015

def year_to_index(year):
    """Transforms a year value to an index in the range given by
       the global variables year_min and year_max"""
    assert(year>=year_min and year <= year_max)  # Error if year is outside range
    return year-year_min

def initialize_temperature_dict(file_name):
    """This code initializes the empty array as values to particular keys in the dictionary """
    temp_new = open(file_name,'r')
    lines = temp_new.readlines()
    temp_dict = {}
    for line in lines:
        station_id = line[:11]# Takes the 11 characters of the station id present in the first column
        empty_array = np.empty(shape=(116,12))#Creates an empty array of the particular size
        empty_array[:] = np.NaN #Fills the array with Not-a-number
        temp_dict[station_id] = empty_array
        temp_new.close()
    return temp_dict




def read_temperature_dict(dat_file,dict):
    import re
    input = open(dat_file, 'r')
    temp_read_data = input.readlines()
    input.close()
    pos_temp = range(19,115,8) #Because the temp data start from pos 19 and end at pos 115. A new temp is seen every 8 pos.
    pattern_temp_miss = re.compile(r'-[9]{4}')
    for line in temp_read_data:
        match_ID = line[0:11]
        match_year = int(line[11:15])
        if match_year >= year_min:
            year_index = year_to_index(match_year)
            month_count = 0
            pos_temp_count = 0
            for position in pos_temp:
                read_temp = line[pos_temp[pos_temp_count]:pos_temp[pos_temp_count]+5]
                if pattern_temp_miss.match(read_temp):
                    pass
                    pos_temp_count += 1
                else:
                    dict[match_ID][year_index,month_count] = float(read_temp)/100.0
                    pos_temp_count += 1
                month_count += 1
        else:
            pass
    return dict

def create_id_to_coordinate_dict(inv_file):
    """Creates the coordinates for the particular weather stations """
    input = open(inv_file, 'r')
    coordinate_data = input.readlines()
    input.close()
    cordi_dict = {}
    for values in coordinate_data:
        x = values.split()
        stat_id =  x[0]
        stat_log = x[1]
        stat_lat = x[2]
        cordi_dict[stat_id] = (stat_log , stat_lat)#Gives keys and assign logitude and latitude to it.
    return cordi_dict



def calculate_temperature_reference(read_temprature_dict, start_year, end_year):
    first_year = year_to_index(start_year)
    last_year = year_to_index(end_year)
    temp_ref_dict = {}
    for i in read_temprature_dict:
        temp_ref_dict[i] = np.nanmean(read_temprature_dict[i][first_year:last_year], axis =0)
    return temp_ref_dict



def calculate_anomalies(read_temperature_dict,calculate_temperature_reference):
    cal_anm_dict = {}
    for keys in read_temperature_dict:
        cal_anm_dict[keys] = np.subtract(read_temperature_dict[keys],calculate_temperature_reference[keys])
    return cal_anm_dict



def anomalies_map_plot(calculate_anomalies, create_id_to_coordinate_dict, start_year, end_year):
    anomalies_dict = {}
    max_temp_anomalies = 0
    min_temp_anomalies = 0

    for station in calculate_anomalies:
        temp_patterns = np.array(calculate_anomalies[station])
        temp_patterns = np.absolute(temp_patterns)
        anomalies_dict.update({station: temp_patterns})
        max_anom = np.amax(temp_patterns[year_to_index(start_year):year_to_index(end_year),:])
        min_anom = np.amin(temp_patterns[year_to_index(start_year):year_to_index(end_year),:])
        if max_anom > max_temp_anomalies:
            max_temp_anomalies = max_anom
        if min_anom > min_temp_anomalies:
            min_temp_anomalies = min_anom

    month_name = {0:"January",1:"February",2:"March",3:"April",4:"May",5:"June",6:"July",7:"August",8:"September",9:"October",10:"November",11:"December"}
    start_index = year_to_index(start_year)
    end_index = year_to_index(end_year)
    plot_count = 0
    year_count = start_year

    for year in range(start_index, end_index):
        for month in range(0,12):
            for station in create_id_to_coordinate_dict:
                latlong = create_id_to_coordinate_dict[station]
                if station in anomalies_dict:
                    lat = latlong[0]
                    long = latlong[1]
                    temp_dev = anomalies_dict[station][year,month]
                    plt.scatter(long,lat,c = temp_dev, vmin = min_temp_anomalies, vmax = max_temp_anomalies, edgecolor = "none")
            plt.axis([-180, 180, -90, 90])
            plt.xlabel('longitude')
            plt.ylabel('latitude')
            plt.title( month_name[month] + '/' + str(year_count))
            plt.savefig('temperature_anomalies_map_%d.png' % plot_count)
            plt.clf()
            plt.show()
            plot_count += 1
            if month == 11:
                year_count += 1
    return anomalies_dict

def anomalies_vs_time_plot(anomalies_dict):
    stations = []
    years = []
    for station in anomalies_dict:
        years = []
        for year in anomalies_dict[station]:
            months = []
            for months in year:
                months.append(months)
            years.append(np.average(months))
        stations.append(years)
    np.reshape(stations,(len(stations),len(years)))

    x = range(len(years))
    y = np.nanmean(stations, axis=0)

    plt.plot(x,y)
    plt.xlabel("years")
    plt.ylabel("average anomaly over months in year for all weather stations")
    plt.title("anomalies over time")
    plt.savefig("anomalies_vs_time.png")
    plt.clf()
    plt.show()



