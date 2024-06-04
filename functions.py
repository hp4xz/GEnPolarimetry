
def Spinup(t, P_0, P_inf, g_sc):
    import numpy as np
    return ((P_0 - P_inf)*np.exp(-1.0*g_sc*t) + P_inf)

def Fit(test,time):
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.optimize import curve_fit
    timeFit=np.linspace(time[0],time[len(time)-1],len(time)*200)
    init_vals = [test[0], test[ len(test)-1 ], 3.0] # [P_0, P_inf, g_sc]
    best_vals, covar = curve_fit(Spinup, time, test, p0=init_vals)
    yFit = Spinup(time, best_vals[0], best_vals[1], best_vals[2])
    yFitLine = Spinup(timeFit, best_vals[0], best_vals[1], best_vals[2])
    P_0 = str(np.round(best_vals[0], 3)) + " +/- " + str(np.round(np.sqrt(abs(covar[0][0])), 3)) 
    P_inf = str(np.round(best_vals[1], 3)) + " +/- " + str(np.round(np.sqrt(abs(covar[1][1])), 3)) 
    T_sc = str(np.round(1.0/best_vals[2], 3)) + " +/- " + str(np.round(np.sqrt(abs(covar[2][2]))/( best_vals[2]*best_vals[2],3)) ) 
    print("P_0 = ", P_0, '\n', "P_inf = ", P_inf, '\n',"T_sc = ", T_sc)
    return P_0,P_inf,T_sc,yFitLine,timeFit,time,test, yFit

def MakeTime(data):
    import numpy as np
    #print(data)
    zeroYear=int(str(data[0][0])[:2])
    zeroMonth=int(str(data[0][0])[2:4])
    zeroDay=int(str(data[0][0])[4:6])
    zeroHour=int(str(data[0][0])[6:8])
    zeroMinute=int(str(data[0][0])[8:10])
    relArray=np.empty(0)
    for i in range(0,len(data)):
        #print(i)
        year=int(str(data[i][0])[:2])
        month=int(str(data[i][0])[2:4])
        day=int(str(data[i][0])[4:6])
        hour=int(str(data[i][0])[6:8])
        minute=int(str(data[i][0])[8:10])
        if (month-zeroMonth)%2 == 0:
            relArray= np.append([relArray],(year-zeroYear)*24*365.2422+(month-zeroMonth)*30.5*24+(day-zeroDay)*24+hour-zeroHour+(minute-zeroMinute)/60)
        elif (month-zeroMonth)%2 != 0:
            relArray= np.append([relArray],(year-zeroYear)*24*365.2422+((month-zeroMonth)*30.5+.5)*24+(day-zeroDay)*24+hour-zeroHour+(minute-zeroMinute)/60)

    return (relArray,int(str(data[0][0])[:10]))   

def MakeData(data):
    import numpy as np
    pcup=np.empty(0)
    pcdown=np.empty(0)
    tcup=np.empty(0)
    tcdown=np.empty(0)
    for i in range(0,len(data)):
        pcup=np.append([pcup],data[i][1])
        pcdown=np.append([pcdown],data[i][2])
        tcup=np.append([tcup],data[i][3])
        tcdown=np.append([tcdown],data[i][4])
    TC=(tcup+tcdown)/2
    PC=(pcup+pcdown)/2
    return(PC,TC,pcup,pcdown,tcup,tcdown)
def MakeFullData(data):
    import numpy as np
    pcup=np.empty(0)
    pcdown=np.empty(0)
    tcUPup=np.empty(0)
    tcUPdown=np.empty(0)
    tcDOWNup=np.empty(0)
    tcDOWNdown=np.empty(0)
    for i in range(0,len(data)):
        pcup=np.append([pcup],data[i][1])
        pcdown=np.append([pcdown],data[i][2])
        tcUPup=np.append([tcUPup],data[i][3])
        tcUPdown=np.append([tcUPdown],data[i][4])
        tcDOWNup=np.append([tcDOWNup],data[i][5])
        tcDOWNdown=np.append([tcDOWNdown],data[i][6])
    TCUP=(tcUPup+tcUPdown)/2
    TCDOWN=(tcDOWNup+tcDOWNdown)/2
    PC=(pcup+pcdown)/2
    return(PC,TCDOWN,TCUP,pcup,pcdown,tcDOWNup,tcDOWNdown,tcUPup,tcUPdown)

def Refine(y):
    import numpy as np
    index=np.empty(0)
    x=np.empty(0)
    minimum=30
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

