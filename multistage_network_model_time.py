# -*- coding: utf-8 -*-
"""
Created on Mon Dec 08 15:24:48 2014

@author: Conor D. Cox For paper Pronounced differences in signal processing 
and synaptic plasticity between piriform-hippocampal network stages:
A prominent role for adenosine Figure 9 last panel
We are improving the code and updates will be posted to github.com/cdcox/multistage_network
This is a multi-stage network model to demonstrate the changes a signal goes through
due to regional facilitation in response to theta firing.
For questions please contact Conor D Cox at cdcox1@gmail.com
"""

from brian2 import *
import numpy as np

def smooth(y, box_pts):
    '''used to smooth regions of each graph against random variations'''
    box = np.ones(box_pts)/box_pts
    y_smooth = convolve(y, box, mode='same')
    return y_smooth
    
def string_assemble(data):
    '''takes input graphs, smooths the, then finds a derivative point far enough away
    for derivative to be 'meaningful' then builds piecewise function based off critical points'''
    data=[float(x) for x in data if x!='']
    data=data[:100]
    data=smooth(data,10)
    data=data[:-10]
    derivative=data[7:]-data[:-7]
    edges=np.where((derivative>0)[1:]-(derivative>0)[:-1])
    output_str='0'
    first_edge=0
    for NN,number in enumerate(edges[0]):
        fit=polyfit(np.arange(first_edge,number)*200,data[first_edge:number],1)
        output_str=output_str+r'+(('+str(fit[0])+r')*t/ms+'+str(fit[1])+r')'+r'*(t/ms>'+str(first_edge*200)+r')*(t/ms<='+str(number*200)+r')'
        first_edge=number
    fit=polyfit(np.arange(first_edge,len(data))*200,data[first_edge:len(data)],1)
    output_str=output_str+r'+(('+str(fit[0])+r')*t/ms+'+str(fit[1])+r')'+r'*(t/ms>'+str(first_edge*200)+r')'
    return(output_str)
    
def build_model_string(equation,firing_strength):
    '''Takes piecewise functions and sets them in a format appropriate for Brian2'''
    string_out=r'scaler='+equation+'''
    x=x+1 
    ge+='''+str(firing_strength)+'''*scaler/100.*mV'''
    return(string_out)

'''Opens book of initial graphs and clean them up for use in Brain2'''    
book = np.genfromtxt(r'5hzresponsedata.csv',delimiter=',',skiprows=3)

#############################################################################
LOT_data=np.array(book[0:96,1])
LOT_str=string_assemble(LOT_data)
 #############################################################################   
AP_data=np.array(book[0:96,2])
AP_str=string_assemble(AP_data)
#############################################################################
PE_data=np.array(book[0:96,3])
PE_str=string_assemble(PE_data)
#############################################################################
LPP_data=np.array(book[0:96,4])
LPP_str=string_assemble(LPP_data)
#############################################################################
DG_data=np.array(book[0:96,5])
DG_str=string_assemble(DG_data)
#############################################################################
CA3_recurrent_data=np.array(book[0:96,6])
CA3_str=string_assemble(CA3_recurrent_data)
#############################################################################
CA1_data=np.array(book[0:96,6])
CA1_str=string_assemble(CA1_data)
#############################################################################
'''Define connectivity matrix based on IO curve values'''
names=['LOT','ASSN','PE','LPP','MF','CA1','CA3',]
x=[1.8,1.8,1.8,1.8,-1,.2,.2]
x=[1-(stuff)/3.1/1.7 for stuff in x]
fix=.15
con_list=np.array(x)*fix
'''baseline equation of excitatory and inhibitory decays'''
eqs = '''
dv/dt = ((ge+gi)/(20.*ms)-(v+60.*mV)/(20.*ms)) : volt
dge/dt = -ge/(5.*ms) : volt
dgi/dt = -gi/(10.*ms) : volt
'''

'''Book keeping to set up model'''
generalized_variable_model='''x:1
                              scaler:1'''
strength=3.2 #excitatory strength
'''Setting up olfactory bulb spiking'''
start=0
indices=[]
times=[]

for i in range(200):
    for j in range(200):
        indices.append(j)
        times.append(start)
    start+=200
indices = array(indices)
times = array(times - 0.5*defaultclock.dt/ms)*ms #moving spiking times to fix Brian2 bug
Olfactory_bulb = SpikeGeneratorGroup(200, indices, times)

'''Each region is set up with its number of cells, equation of firing strength change over time
and connectivity and after hyper-polarization'''
Piriform=NeuronGroup(200, eqs, threshold='v>-50*mV', reset='v=-60*mV')
Piriform.v = -60*mV
model_string=build_model_string(LOT_str,strength)
LOT = Synapses(Olfactory_bulb, Piriform, model=generalized_variable_model,pre=model_string)
LOT.connect(True, p=con_list[0])
##############################################################################
model_string=build_model_string(AP_str,strength)
CP_P=Synapses(Piriform, Piriform, model=generalized_variable_model,pre=model_string)                                            
CP_P.connect(True,  p=con_list[1])
CP_P.delay=1*ms  
##########################################################################################
CP_ah=Synapses(Piriform,Piriform,pre=r'gi-=4*mV')
for i in range(200):
    CP_ah.connect(i,i)
#############################################################################################
Entorhinal=NeuronGroup(200, eqs, threshold='v>-50*mV', reset='v=-60*mV') 
Entorhinal.v = -60*mV
model_string=build_model_string(PE_str,strength)
CP_e= Synapses(Piriform, Entorhinal,model=generalized_variable_model,pre=model_string)
CP_e.connect(True, p=con_list[2])
#############################################################################
CE_ah=Synapses(Entorhinal,Entorhinal,pre=r'gi-=4*mV')
for i in range(200):
    CE_ah.connect(i,i)
############################################################################
Dentate=NeuronGroup(200, eqs, threshold='v>-50*mV', reset='v=-60*mV')
Dentate.v = -60*mV
model_string=build_model_string(LPP_str,strength)
LPP = Synapses(Entorhinal, Dentate,model=generalized_variable_model,pre=model_string)
LPP.connect(True, p=con_list[3])
##############################################################################
CD_ah=Synapses(Dentate,Dentate,pre=r'gi-=4*mV')
for i in range(200):
    CD_ah.connect(i,i)
############################################################################
CA3=NeuronGroup(400, eqs, threshold='v>-50*mV', reset='v=-60*mV')
CA3.v = -60*mV
model_string=build_model_string(DG_str,10)
MF = Synapses(Dentate, CA3,model=generalized_variable_model,pre=model_string)
MF.connect(True, p=con_list[4])
###########################################################################
'''simulation of simplified CA3 inhibition'''
C3_ah=Synapses(CA3,CA3,pre=r'gi-=4*mV')
for i in range(400):
    C3_ah.connect(i,i)
    
C3_recurrents=Synapses(CA3,CA3,pre=r'gi-=2*mV')
C3_recurrents.connect(True,p=.15)
C3_recurrents.delay=2*ms
######################################################################
model_string=build_model_string(CA3_str,strength)
C3_c3 = Synapses(CA3, CA3,model=generalized_variable_model,pre=model_string)
C3_c3.connect(True,p=con_list[5])
C3_c3.delay=1*ms                              
############################################################################
CA1=NeuronGroup(200, eqs, threshold='v>-50*mV', reset='v=-60*mV')
CA1.v = -60*mV
model_string=build_model_string(CA1_str,strength)
SC = Synapses(CA3, CA1, model=generalized_variable_model,pre=model_string)
SC.connect(True, p=con_list[6])
###########################################################################
C1_ah=Synapses(CA1,CA1,pre=r'gi-=4*mV')
for i in range(200):
    C1_ah.connect(i,i)

'''clear up any open figures'''
cla()   # Clear axis
clf()   # Clear figure
close() # Close a figure window
'''inform user that simulation is under way'''
print 'initialized'
'''set up monitoring for each region'''
Mp = PopulationRateMonitor(Piriform)
Me = PopulationRateMonitor(Entorhinal)
Md = PopulationRateMonitor(Dentate)
Mc = PopulationRateMonitor(CA3)
Mc1 = PopulationRateMonitor(CA1)

'''example state-monitors showing to monitor values of individual cells and 
values of synapse groups see Brian documentation for more details about how to use these'''
cell_monitor= StateMonitor(CA3,'v',record=CA3)  
synapse_monitor= StateMonitor(C3_c3,'scaler',record=C3_c3[0,:])  
'''run the model'''
run(20*second)

'''average bins into firing, normalize and plot'''
rate_sump=[]
rate_sume=[]
rate_sumd=[]
rate_sumc=[]
rate_sumc1=[]
for zed in range(len(Mp.rate)/2000):
    rate_sump.append(np.sum(Mp.rate[((zed-1)*2000+5):(zed*2000+5)]/Hz)) 
    rate_sume.append(np.sum(Me.rate[((zed-1)*2000+5):(zed*2000+5)]/Hz)) 
    rate_sumd.append(np.sum(Md.rate[((zed-1)*2000+5):(zed*2000+5)]/Hz)) 
    rate_sumc.append(np.sum(Mc.rate[((zed-1)*2000+5):(zed*2000+5)]/Hz)) 
    rate_sumc1.append(np.sum(Mc1.rate[((zed-1)*2000+5):(zed*2000+5)]/Hz))
rate_sump=rate_sump/np.max(rate_sump)
rate_sume=rate_sume/np.max(rate_sume)
rate_sumd=rate_sumd/np.max(rate_sumd)
rate_sumc=rate_sumc/np.max(rate_sumc)
rate_sumc1=rate_sumc1/np.max(rate_sumc1)
plot(Mp.t[::2000]/ms, rate_sumc1,'r',Mp.t[::2000]/ms,rate_sumc,'g',Mp.t[::2000]/ms,rate_sump,'b',Mp.t[::2000]/ms,rate_sumd,'k',Mp.t[::2000]/ms,rate_sume,'c')
show()