
import os
import sys
sys.path.insert(0, "..")

##tutaj biblioteki do liczenia i wyświetlania
import numpy as np
from scipy.linalg import hankel

from scipy.integrate import odeint
from scipy import linalg
from scipy.signal import lfilter
import scipy.io as sio

##do stuktury p -> dataclass
from dataclasses import dataclass
from typing import List

# do odbsługi serwera
from opcua import Client
from opcua import ua
# do cyklu działania programu
import sched, time
s = sched.scheduler(time.time, time.sleep)

Step_response = sio.loadmat(os.path.dirname(os.path.abspath(__file__))+"\step_response_1sec_final_4ord_197.mat")

print(Step_response.keys())

@dataclass
class D_struct:
    sr: List[np.double]
    p: int
    m: int
    la: int
    y: float
    v: List[np.double]
    r: float
    F: List[np.double]
    G: List[np.double]
    k: List[np.double]
    u: np.double

def dmc(D):
    #print("DMC")

    #size of step response
    #ilosc danych w step response ukladu
    N=np.size(D.sr)

    # predition horizon parameter
    # parametr horyzontu przewidywania
    Dc=D.p
    #print(Dc)

    # first run definition
    # definicje pierwszego uruchomienia
    if not D.v:
        # number of past inputs to keep
        # liczba poprzednich wejść
        n=N-Dc
        # storage for past input
        # wektor poprzednich wejść
        D.v=np.zeros(n)
        # matrix to calculate free response from past input
        # macierz do liczenia skladowej swobodnej z poprzednich wejść
        x=D.sr[0:n]
        # tutaj ze wzoru (10) Ku
        D.F=hankel(D.sr[1:Dc+1],D.sr[Dc:N])-np.tile(x,[Dc,1])
        
        # dynamic matrix
        # macierz dynamiczna M
        first_row = D.sr[1]*np.eye(1,D.m)
        D.M=linalg.toeplitz(D.sr[0:Dc],first_row)
        
        # calculate DMC gain
        # liczenie wzmocnienia DMC
        #print(D.la*np.eye(D.m))
        R=np.linalg.cholesky(np.dot(np.transpose(D.M),D.M) + D.la*np.eye(D.m))
        # print(R) z jakiegoś powodu wychodzi tranponowane R
        R=np.transpose(R)
        rightRnaDG = np.linalg.lstsq(np.transpose(R),np.transpose(D.M),rcond=None)
        #macierz współczynników K jest obliczana raz w fazie obliczania regulatora
        K=np.linalg.lstsq(R,rightRnaDG[0],rcond=None)
        
        # only the first input will be used
        # tylko pierwsze używane
        #D.k=K[1,:]  
        #tutaj ze wzoru (10) Ke
        D.k=K[0][0] 
        D.u=0
    #end of first run definition
    #koniec definicji pierwszego uruchomienia

    # free response
    # skadowa swobodna
    f=np.dot(D.F,D.v)
    
    # reference sp
    #referencja wartości zadanej
    nr=np.size(D.r)
    if(nr>=Dc):
        ref=D.r[0:Dc]
    else:
        ref=[D.r]
        ref.extend(D.r[-1]+np.zeros(Dc-nr))
    

    ref[0]=D.y

    # DMC input change
    # obliczenie nowego du z 
    du=np.dot(D.k,(ref-f-D.y))

    # print("ref")
    # print(ref)
    # print("f")
    # print(f)
    # print("u")
    # print(u)

    # past input change for next step
    # przesunięcie wektora poprzednich wejść
    prevDvWOlast = D.v[:-1]
    prevDvWOlast=np.append(du,prevDvWOlast)
    D.v = prevDvWOlast.tolist()
    
    # next input
    # następne u
    D.u=D.u+du

    # so that the output is in <0,1>
    # ograniczenie parametrów do przedziału <0,1>
    # anti-wind-up
    if D.u > 1 :
        D.u = 1
    if D.u < 0 :
        D.u = 0

    return D

def dmc_cycle(sc,client,D):
    
    print("")

    # print("D.p")
    # print(D.p)
    # print("D.m")
    # print(D.m)
    # print("D.la")
    # print(D.la)

    var_sp = client.get_node('ns=3;s="DB_Instance_of_OPC_UA_DB"."OPC Data"."ua_setpoint"')
    #var_get = var_sp.get_value()
    
    # wstawienie aktualnego setpoint 
    D.r=[] 
    D.r.append([var_sp.get_value()]) 

    print("setpoint - reference (set point)")
    print(D.r)
    
    var_DMC_submit_parameters = client.get_node('ns=3;s="DB_Instance_of_OPC_UA_DB"."OPC Data"."ua_DMC_submit_parameters"')

    #podczas zmiany parametrów regulacji
    if var_DMC_submit_parameters.get_value() == True:
        D.v=[]
        
        param_p = client.get_node('ns=3;s="DB_Instance_of_OPC_UA_DB"."OPC Data"."ua_p"')
        param_m = client.get_node('ns=3;s="DB_Instance_of_OPC_UA_DB"."OPC Data"."ua_m"')
        param_la = client.get_node('ns=3;s="DB_Instance_of_OPC_UA_DB"."OPC Data"."ua_la"')
        
        D.p = param_p.get_value()
        D.m = param_m.get_value()
        D.la = param_la.get_value()
    
    print("D.p")
    print(D.p)
    print("D.m")
    print(D.m)
    print("D.la")
    print(D.la)


    # #jeśli dmc zostaje uruchomiony
    # var_DMC_change = client.get_node('ns=3;s="DB_Instance_of_OPC_UA_DB"."OPC Data"."ua_DMC_change"')
    # print("DMC button pressed")
    # print(var_DMC_change.get_value())
    # if(var_DMC_change.get_value() == True):
    #     D.u = 0
    # print(D.u)

    D = dmc(D)

    print("D.u - current input")
    print(D.u)
    
    var_u = client.get_node('ns=3;s="DB_Instance_of_OPC_UA_DB"."OPC Data"."ua_u"')
    dv = ua.DataValue(ua.Variant(D.u, ua.VariantType.Float))
    var_u.set_value(dv)
    
    var_y = client.get_node('ns=3;s="DB_Instance_of_OPC_UA_DB"."OPC Data"."ua_y"')
    # oddanie rzeczywistego wyjscia
    D.y  = var_y.get_value()
    
    print("D.y - current mrasurement")
    print(D.y)
    

    #s.enter(0.1, 1, dmc_cycle, (sc,client,D))
    s.enter(1, 1, dmc_cycle, (sc,client,D))

#do uruchomienia
if __name__ == "__main__":
    print("MAIN -> uruchamiam regulator DMC")
    D=D_struct(
        np.transpose(Step_response['y'])[0] # D.sr - unit step response data
        ,43 #56 #20 # D.p  - prediction horizon 20 to 70
        ,101 #133 #5 # D.m  - dynamic (moving) horizon
        ,28.38#21.5#43#28.38#14.19 #10 # D.la - performance criterion weight
        ,0  # D.y  - current mrasurement
        ,[] # D.v  - past input, initially empty
        ,[]  # D.r  - reference (set point)
        ,[] # D.F  - matrix to calculate free response, set by the initial call
        ,[] # D.M  - dynamic matrix, set by the initial call
        ,[] # D.k  - DMC gain, set by the initial call
        ,None # D.u  - current input, initially 0
        )
 
    print("len(D.sr)")
    print(len(D.sr))
    print(D.sr)

    client = Client("opc.tcp://192.168.1.5:4840")
    try:
        client.connect()

        # Client has a few methods to get proxy to UA nodes that should always be in address space such as Root or Objects
        root = client.get_root_node()
        print("Objects node is: ", root)

        # Node objects have methods to read and write node attributes as well as browse or populate address space
        print("Children of root are: ", root.get_children())
        
        s.enter(0, 1, dmc_cycle,  (s,client,D,))
        s.run()

    finally:
        client.disconnect()

