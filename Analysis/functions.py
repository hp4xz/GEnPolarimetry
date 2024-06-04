def Calibrate(afpDate,firstNMRDate,secondNMRDate,percent_per_khz,epicsData):
    import re
    import numpy as np
    
    trash, fShift, bShift=FindEPRShifts(afpDate)
    fit=np.load('NMRDataFits/FitsInRange.npy')
   
    firstNMR=np.where(fit[0].astype(int)==firstNMRDate)[0][0]
    secondNMR=np.where(fit[0].astype(int)==secondNMRDate)[0][0]
    #Get percent by multiplying front and back shift by %/khz
    fPer=fShift*percent_per_khz   
    bPer=bShift*percent_per_khz 
    
    f_pc_const=fPer/fit[1][firstNMR]
    b_pc_const=bPer/fit[1][secondNMR]
    
    f_us_const=fPer/fit[2][firstNMR]
    b_us_const=bPer/fit[2][secondNMR]
    
    f_ds_const=fPer/fit[3][firstNMR]
    b_ds_const=bPer/fit[3][secondNMR]
    
    pc_const=np.round((f_pc_const+b_pc_const)/2,4)
    us_const=np.round((f_us_const+b_us_const)/2,4)
    ds_const=np.round((f_ds_const+b_ds_const)/2,4)
    print(' Format of output is: EPR date, pc_const, us_const, ds_const')
    if epicsData!="none":
        epicsFirst1=str(firstNMRDate)[:8]
        epicsFirst2=str(firstNMRDate)[8:]
        epicsFirst=epicsFirst1+'_'+epicsFirst2

        #print(epicsFirst)
        epicsSecond1=str(secondNMRDate)[:8]
        epicsSecond2=str(secondNMRDate)[8:]
        epicsSecond=epicsSecond1+'_'+epicsSecond2

        firstTemp=np.empty(0)
        secondTemp=np.empty(0)
        for i in epicsData:   
            firstTemp=np.append(firstTemp,np.round(GrabAround(epicsFirst,i),4))
            secondTemp=np.append(secondTemp,np.round(GrabAround(epicsSecond,i),4))
    
        return int(re.search(r'\d+', afpDate).group()), pc_const,us_const,ds_const,firstTemp,secondTemp
    else:
         return int(re.search(r'\d+', afpDate).group()), pc_const,us_const,ds_const,0,0
    
def DatesInRange(start,end):
    from datetime import datetime, timedelta
    import os
    import numpy as np
    NMRFolders=os.listdir('NMRData')
    foldersRefined=np.empty(0)
    inRange=np.empty(0)
    inRangeIndex=np.empty(0)
    for i in range(0,len(NMRFolders)):
        elements=NMRFolders[i].split('_')
        combined=elements[2]+'_'+elements[3]
        foldersRefined=np.append(foldersRefined,combined)
    start_date = datetime.strptime(start, "%Y%m%d_%H%M%S")
    end_date = datetime.strptime(end, "%Y%m%d_%H%M%S")
    for i in range(0,len(foldersRefined)):
        date = datetime.strptime(foldersRefined[i], "%Y%m%d_%H%M%S")
        if start_date <= date <= end_date:
            inRangeIndex=np.append(inRangeIndex,i)
            inRange=np.append(inRange,foldersRefined[i])
    return inRangeIndex.astype(int),inRange
def FindEPRShifts(afpDate):
    import numpy as np
    file=np.loadtxt(afpDate)
    data=np.transpose(file)
    first=list(ordered_cluster(data[1], 10))
    second=Refine(first)
    return second
def FindNMRsB4andAft(date, date_array):
    from datetime import datetime
    import numpy as np
    import os
    date = datetime.strptime(date, "%Y%m%d%H%M")
    before_index = None
    after_index = None


    for i, d in enumerate(date_array):
        sameMin=0
        found=0

        current_date = datetime.strptime(d, "%Y%m%d_%H%M%S")

        
        rounded_date = current_date.replace(second=0)
        current_date=rounded_date
        if current_date < date and current_date!= date:
          
            before_index = i
            found=1
         
        elif current_date==date:
            before_index=i
            after_index=i+1
           
            break
       
        elif current_date > date:
            
            after_index = i
            
            break

    return before_index, after_index
def FitAll(indicesInRange,dates):
    print('Format of output: dates,pc,us,ds,pcusx,pcdsx,pcusy,pcdsy,dsusx,dsdsx,dsusy,dsdsy,ususx,usdsx,ususy,usdsy,FWHM 1 or 0')
    import numpy as np
    import os
    from math import sqrt
    fitData=np.empty(0)
    cov=[]
    for i in range(0,len(indicesInRange)):
        try:
            #print('1')
            fwhm='true'
            time=dates[i].replace('_','')
            oneSweep=GrabSweepAtIndex(indicesInRange[i])
            #print('2')
            dsusx,params,cov_dsusx = FitLorentzian(0,1,oneSweep,i)
            #if params[3]>
            dsdsx,params,cov_dsdsx = FitLorentzian(0,2,oneSweep,i)
            dsusy,params,cov_dsusy = FitLorentzian(0,3,oneSweep,i)
            dsdsy,params,cov_dsdsy = FitLorentzian(0,4,oneSweep,i)

            ususx,params,cov_ususx = FitLorentzian(1,1,oneSweep,i)
            usdsx,params,cov_usdsx = FitLorentzian(1,2,oneSweep,i)
            ususy,params,cov_ususy = FitLorentzian(1,3,oneSweep,i)
            usdsy,params,cov_usdsy = FitLorentzian(1,4,oneSweep,i)

            pcusx,params,cov_pcusx = FitLorentzian(2,1,oneSweep,i)
            pcdsx,params,cov_pcdsx = FitLorentzian(2,2,oneSweep,i)
            pcusy,params,cov_pcusy = FitLorentzian(2,3,oneSweep,i)
            pcdsy,params,cov_pcdsy = FitLorentzian(2,4,oneSweep,i)

            #print('3')

            pcust = sqrt(pcusx**2+pcusy**2)
            pcdst = sqrt(pcdsx**2+pcdsy**2)

            usust = sqrt(ususx**2+ususy**2)
            usdst = sqrt(usdsx**2+usdsy**2)

            dsust = sqrt(dsusx**2+dsusy**2)
            dsdst = sqrt(dsdsx**2+dsdsy**2)
            cov.append([cov_pcusx,cov_pcdsx,cov_pcusy,cov_pcdsy,cov_dsusx,cov_dsdsx,cov_dsusy,cov_dsdsy,cov_ususx,cov_usdsx,cov_ususy,cov_usdsy])


            pc=(pcust+pcdst)/2

            us=(usust+usdst)/2

            ds=(dsust+dsdst)/2

            fitData= np.append(fitData,[int(time),pc,us,ds,pcusx,pcdsx,pcusy,pcdsy,dsusx,dsdsx,dsusy,dsdsy,ususx,usdsx,ususy,usdsy],axis=0)
        
        except Exception as e:
            print('Failed to Fit Dataset: ',i,' ', str(e))
    return fitData,cov
def FitLorentzian(n,q,oneSweep,i):
    error=0
    import numpy as np
    from scipy.optimize import curve_fit
    
    # Define the Lorentzian function
    def lorentzian(x, ampl, center, fwhm):
        return ampl / (1 + ((x - center) / (0.5 * fwhm))**2)
    x_data = oneSweep[n][0]
    y_data = oneSweep[n][q]
    COV=[]
    
    ##added this to help with wonky fits during fringe and chicago running
    ymean=np.mean(y_data)
    if(abs(ymean)>20):
        y_data-=ymean
        
        
    ymax=np.amax(y_data)
    ymaxarg=np.argmax(y_data)
    ymin=np.amin(y_data)
    yminarg=np.argmin(y_data)

    if abs(ymax)>abs(ymin):

        initial_params=[ymax,x_data[ymaxarg],1]
    else:
        initial_params=[ymin,x_data[yminarg],1]



    # Fit the Lorentzian function to the data
    try:
        
        params, COV = curve_fit(lorentzian, x_data, y_data,p0=initial_params)

        # Extract the height from the fit parameters
        height = params[0]
        
    except:
        params=[0,0,0]
        height=0
       
        error=1
    # Print the height
    #print("Height =", height)
    return(height,params,COV)

def GrabAroundGPT(date, data):
    import time
    import numpy as np
    from datetime import timedelta
    use_date = time.strptime(date, "%Y%m%d_%H%M%S")
    start = use_date - timedelta(minutes=2)
    end = use_date + timedelta(minutes=2)
    
    dates = data[0]
    ydat = data[1].astype(float)
    
    mask = np.isin(dates, [start, end])
    in_range_indices = np.where(mask)[0]
    refined = ydat[in_range_indices]
    
    return refined

def GrabAround(date,data):
    from datetime import datetime, timedelta
    import os
    import numpy as np
    refined=np.empty(0)
    inRange=np.empty(0)
    inRangeIndex=np.empty(0)

    useDate=datetime.strptime(date, "%Y%m%d_%H%M%S")
    start = useDate - timedelta(minutes=2)
    end = useDate + timedelta(minutes=2)
    yep=np.empty(0)
    h=0
    ydat=data[1].astype(float)
    #print(ydat[0])
    for i in range(0,len(data[1])):
        thisDate = datetime.strptime(data[0][i], "%Y%m%d_%H%M%S")
        
        if start <= thisDate <= end:
            #print('found ',i)
            yep=np.append(yep,ydat[i])
    #print(yep)
            
        
    
    return(np.mean(yep))  
def GrabSettingsAtIndex(j):
    import os
    import numpy as np
    import pandas as pd
    import re
    NMRFolders=os.listdir('NMRData')
    j=30
    name='NMRData/'+NMRFolders[j]+'/Downstream Coil/Settings.dat'
    with open(name, 'r') as f:
        # read the contents of the file into a list
        lines = f.readlines()

    # remove newline characters from each element
    lines = [line.strip() for line in lines]

    # print the list

    result=[]
    for string in lines[4:11]:
        # split the string into key and value
        key, value = string.split(':')
        # strip leading and trailing whitespace from key and value
        key = key.strip()
        value = value.strip()
        # convert value to integer
        value = float(value)
        # add key and value to result list
        result.append([key, value])
    return result
def GrabSweepAtIndex(j):
    import os
    import numpy as np
    NMRFolders=os.listdir('NMRData')

    dates=np.empty(0)
    allSweeps=np.empty((0,0,0))

    downstreamFolder=os.listdir('NMRData/'+NMRFolders[j]+'/Downstream Coil')
    upstreamFolder=os.listdir('NMRData/'+NMRFolders[j]+'/Upstream Coil')
    pumpingchamberFolder=os.listdir('NMRData/'+NMRFolders[j]+'/Pumping Chamber')


    
    for i, s in enumerate(downstreamFolder):
        if s.startswith('0001') and s.endswith('.dat'):
            DSindex = i
            break

    for i, s in enumerate(upstreamFolder):
        if s.startswith('0001') and s.endswith('.dat'):
            USindex = i
            break

    for i, s in enumerate(pumpingchamberFolder):
        if s.startswith('0001') and s.endswith('.dat'):
            PCindex = i
            break

    singleDSData=np.expand_dims(np.transpose(np.loadtxt('NMRData/'+NMRFolders[j]+'/Downstream Coil/'+downstreamFolder[DSindex])),axis=0)
    singleUSData=np.expand_dims(np.transpose(np.loadtxt('NMRData/'+NMRFolders[j]+'/Upstream Coil/'+upstreamFolder[USindex])),axis=0)
    singlePCData=np.expand_dims(np.transpose(np.loadtxt('NMRData/'+NMRFolders[j]+'/Pumping Chamber/'+pumpingchamberFolder[PCindex])),axis=0)
    dates=np.append(dates,NMRFolders[j])

    oneSweep=np.concatenate((singleDSData,singleUSData,singlePCData),axis=0)
    #the first index of oneSweep is the downstream data, second is the upstream data, third is the pumping chamber data
    return oneSweep






def LoadEpics():
    import numpy as np
    import sys
    import os
    sys.path.append(sys.path[0]+'/../../')
    tc1=np.load('Epics_NP_Arrays/tc1.npy')
    tc2=np.load('Epics_NP_Arrays/tc2.npy')
    tc3=np.load('Epics_NP_Arrays/tc3.npy')
    tc4=np.load('Epics_NP_Arrays/tc4.npy')

    pc1=np.load('Epics_NP_Arrays/pc1.npy')
    pc2=np.load('Epics_NP_Arrays/pc2.npy')

    tt1=np.load('Epics_NP_Arrays/tt1.npy')
    tt2=np.load('Epics_NP_Arrays/tt2.npy')
    print('Format: tc1,tc2,tc3,tc4,pc1,pc2,tt1,tt2')
    return tc1,tc2,tc3,tc4,pc1,pc2,tt1,tt2
def Lorentzian(x, ampl, center, fwhm):
        return ampl / (1 + ((x - center) / (0.5 * fwhm))**2)
def ordered_cluster(data, max_diff):
    from statistics import mean
    current_group = ()
    for item in data:
        test_group = current_group + (item, )
        test_group_mean = mean(test_group)
        if all((abs(test_group_mean - test_item) < max_diff for test_item in test_group)):
            current_group = test_group
        else:
            yield current_group
            current_group = (item, )
    if current_group:
        yield current_group
def Refine(y):
    import numpy as np
    index=np.empty(0)
    x=np.empty(0)
    minimum=25
    for i in range(0,len(y)):
        add=np.arange(0,len(y[i]),1)
        if len(add)> minimum:
            x=np.append(x,np.mean(y[i]))
            
        else:
            index=np.append(index,i)
     
    index=index.astype(int)
    ynew=np.delete(y,index)
    
    n=10
    if (len(ynew) == 5):
        for i in range(0,len(ynew)):
            ynew[i]=ynew[i][n:]
            ynew[i]=ynew[i][:-n]
        Front=np.mean(ynew[0][-5:])-np.mean(ynew[1][5:])
        Back=np.mean(ynew[4][5:])-np.mean(ynew[3][-5:])
        print(' Front shift is: ', np.round(Front,2), 'kHz \n', 'Back shift is: ', np.round(Back,2), 'kHz')
    elif (len(ynew)==3):
        for i in range(0,len(ynew)):
            ynew[i]=ynew[i][n:]
            ynew[i]=ynew[i][:-n]
        Front=np.mean(ynew[0][-5:])-np.mean(ynew[1][5:])
        Back=np.mean(ynew[2][5:])-np.mean(ynew[1][-5:])
        print(' Front shift is: ', np.round(Front,2), 'kHz \n', 'Back shift is: ', np.round(Back,2), 'kHz')
    else:
        print('Incorrect Grouping')
        Front=0
        Back=0
    
        
    return(ynew,Front,Back)


