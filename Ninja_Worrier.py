import numpy as np
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import sys

print("************************************************************")
print("**************** ||     MAJOR PROJECT     ||****************")
print("************************************************************")

print('By: Anuj Patel, Prachi Pardhi, Sakshi Pathak')

seismic_zones = {
    "Delhi": "Zone IV",
    "Mumbai": "Zone III",
    "Chennai": "Zone III",
    "Kolkata": "Zone III",
    "Bengaluru": "Zone III",
    "Hyderabad": "Zone III",
    "Pune": "Zone III",
    "Ahmedabad": "Zone III",
    "Jaipur": "Zone II",
    "Surat": "Zone III",
    "Lucknow": "Zone II",
    "Kanpur": "Zone II",
    "Nagpur": "Zone II",
    "Visakhapatnam": "Zone III",
    "Bhopal": "Zone II",
    "Patna": "Zone II",
    "Ludhiana": "Zone II",
    "Agra": "Zone II",
    "Nashik": "Zone II",
    "Vadodara": "Zone III",
    "Faridabad": "Zone II",
    "Madurai": "Zone II",
    "Jamshedpur": "Zone II",
    "Nasik": "Zone II",
    "Coimbatore": "Zone II",
    "Jammu": "Zone IV",
    "Amritsar": "Zone IV",
    "Jodhpur": "Zone IV",
    "Raipur": "Zone III",
    "Kota": "Zone IV",
    "Guwahati": "Zone V",
    "Indore":"Zone II"
}

# =============================================================================
# Input Part
# =============================================================================

height=float(input("Please Enter the Height of each Storey in meters : "))
storey=int(input("Please Enter the Number of Storey : "))
totalHeight = storey*height

bay_X=int(input("Please Enter the Number of Bay in X direction : "))
bay = bay_X
print('')
print('---||---')

widths_bay=[]
print("please input the Widths of bay :-")
lenX=0
for i in range(bay_X):
    print('Bay ',i+1,end=' ')
    a=float(input("Width = "))
    widths_bay.append(a)
    lenX+=a
    
#################################################################################
    
# inputs
# speak('Enter the height of building in meters')
h = totalHeight #height of building
bay_Y = int(input("Enter the number of bays in Y direction: "))
bay_width = {}
lenY=0
for i in range(bay_Y):
    print('Bay ',i+1,end=' ')
    a=float(input("Width = "))
    bay_width[i+1] = a
    lenY+=a

# To find weight
l_c,w_c = input('Enter length and width of column in meters (space separated):').split()
l_c = float(l_c)
w_c = float(w_c)
d_b,w_b = input('Enter depth and width of Beams in meters (space separated):').split()
d_b = float(d_b)
w_b = float(w_b)
slab = float(input('Enter thickness of slab in meters: '))
imposedLoad = float(input('Enter Imposed uniformity distributed floor load (in KN/sq.M): '))




print()
# speak('Enter location of building (City or state)')
location = input("Enter location of building (City or state):  ").title()

print()
building_type = input("Enter type of building (Residential, Commercial, Industrial): ").lower()

print()
print("For your reference:")
print('Type-I : Rock or hard soils ')
print('Type-II : Medium or stiff soils')
print('Type-III : Soft Soils')
soilType = input("Enter the soil type = ").upper()


print()
print("For your reference:")
print('Type-I : RC-MRF building ')
print('Type-II : Rc-Steel composite MRF building')
print('Type-III : MRF building')
frame = input("Enter the type of frame = ").upper()

print()
frameType = input('Enter type of Lateral load resisting system (OMRF/SMRF/OBF/SBF): ').upper()

if frameType == 'OMRF':
    response_reduction_factor = 3.0

if frameType == 'SMRF':
    response_reduction_factor = 5.0
    
if frameType == 'OBF':
    response_reduction_factor = 4.0
    
if frameType == 'SBF':
    response_reduction_factor = 5.0



if location in ["Indore","Delhi", "Mumbai", "Chennai", "Kolkata", "Bengaluru", "Hyderabad", "Pune", "Ahmedabad", "Jaipur", "Surat", "Lucknow", "Kanpur", "Nagpur", "Visakhapatnam", "Bhopal", "Patna", "Ludhiana", "Agra", "Nashik", "Vadodara", "Faridabad", "Madurai", "Jamshedpur", "Nasik", "Coimbatore", "Jammu", "Amritsar", "Jodhpur", "Raipur", "Kota", "Guwahati"]:
    if location in ["Delhi", "Mumbai", "Chennai", "Kolkata", "Bengaluru", "Hyderabad", "Pune"]:
        seismic_zone_factor = 0.16
    elif location in ["Ahmedabad", "Jaipur", "Surat", "Lucknow", "Kanpur", "Nagpur", "Visakhapatnam"]:
        seismic_zone_factor = 0.12
    elif location in ["Bhopal","Indore", "Patna", "Ludhiana", "Agra", "Nashik", "Vadodara"]:
        seismic_zone_factor = 0.10
    elif location in ["Faridabad", "Madurai", "Jamshedpur", "Nasik", "Coimbatore"]:
        seismic_zone_factor = 0.08
    elif location in ["Jammu", "Amritsar", "Jodhpur", "Raipur", "Kota"]:
        seismic_zone_factor = 0.06
    elif location == "Guwahati":
        seismic_zone_factor = 0.05

    if building_type == "residential":
        importance_factor = 1
    elif building_type == "commercial":
        importance_factor = 1.5
    elif building_type == "industrial":
        importance_factor = 1.5
    else:
        print("Invalid building type entered")
        sys.exit()
else:
    print("Location not found")
    sys.exit()
    
R = response_reduction_factor
I = importance_factor
    



# Finding the value of fundamental time period (T)
if frame == 'I':
    T = 0.075 * h**(0.75)
    T=float("{:.2f}".format(T))
if frame == 'II':
    T = 0.080 * h**(0.75)
    T=float("{:.2f}".format(T))
if frame =='III':
    T = 0.085 * h**(0.75)
    T=float("{:.2f}".format(T))
    


# Finding value of (Sa/g) let say it is S
if soilType == 'Rocky' or soilType=='Hard' or soilType == 'I':
    if T<0.40 and T>0:
        S=2.5
    elif T>0.40 and T<4.00:
        S=(1/T)
    elif T>4.00:
        S=0.25
        
if soilType=='Medium Stiff' or soilType == 'II':
    if T<0.55 and T>0:
        S=2.5
    elif T>=0.55 and T<4.00:
        S=(1.36/T)
    elif T>4.00:
        S=0.34
        
if soilType=='Soft' or soilType == 'III':
    if T<0.67 and T>0:
        S=2.5
    elif T>0.67 and T<4.00:
        S=(1.67/T)
    elif T>4.00:
        S=0.42
        


# Finding the Seismic Zone factor (Z)
zone = seismic_zones[location]
if zone == 'Zone II':
    Z = 0.10
if zone == 'Zone III':
    Z = 0.16
if zone == 'Zone IV':
    Z = 0.24
if zone == 'Zone V':
    Z = 0.36
    
    
#### Weight Calculation
numberOfBeams = (bay_X*(bay_Y + 1)) + (bay_Y*(bay_X + 1))
numberOfcolumns = (bay_X + 1)*(bay_Y + 1)
areaOfFloor = (lenY)*(lenX)

weightOfBeams = numberOfBeams*25*d_b*w_b*widths_bay[0]
weightOfSlab = areaOfFloor*slab*25
weightOfcolumns = numberOfcolumns*height*l_c*w_c*25

if imposedLoad<=3:
    liveLoad = (25/100)*imposedLoad*areaOfFloor
elif imposedLoad>3:
    liveLoad = (50/100)*imposedLoad*areaOfFloor
    
weightOfRoof = (weightOfcolumns/2) + weightOfSlab + weightOfBeams + ((25/100)*liveLoad)

######## Total weight of building
totalWeight = ((storey-1)*(weightOfSlab + weightOfcolumns + weightOfBeams + liveLoad)) + weightOfRoof



    

print()
print('>>>>>>>>>> Following are the outputs :- <<<<<<<<<<')
print('Seismic zone factor (Z) = ',Z) # this line prints the value of Z
print('Fundamental time period (T) = ',T) # this line prints the value of T
print('Design acceleration coefficient (Sa/g) = ',S)  # this line prints the value of S
print('Response reduction factor (R) = ',R)  # this line prints the value of R and I
print('Importance factor (I) =',I)

# Formula used to find design horizontal seismic coefficient(Ah)
# here S=Sa/g
A_h = ((Z/2)*S)/(R/I)

# Base Shear
Vb = A_h * totalWeight

print('Design horizontal seismic coefficient (Ah) =',A_h)
print('Design Base Shear Vb =',Vb)


force=[]

hh = 0
sumOffloorWeight = 0
for i in range(1, storey+1):
    hh += height
    hSquare = hh*hh
    print('height of floor',i,'is',hh)
    if i!= storey:
        wj = (weightOfSlab + weightOfcolumns + weightOfBeams + liveLoad)
    else:
        wj =  weightOfRoof
    sumOffloorWeight += wj*hSquare
   
print()
hh2 = 0
for i in range(1, storey+1):
    hh2 += height
    h2Square = hh2*hh2
    print('height of floor',i,'is',hh2)
    if i!= storey:
        wi = (weightOfSlab + weightOfcolumns + weightOfBeams + liveLoad)
    else:
        wi =  weightOfRoof
    Q = ((wi*h2Square)/sumOffloorWeight)*Vb
    Q=float("{:.2f}".format(Q))
    force.append(Q)

print()
print('Forces :',force)
    




##########################################################################################

##### 
    
    
    
    
    
    
    
    
# print('')
# print('---||---')
# force=[]
# print("please input the Forces on Storey :-")
# for i in range(storey):
#     print('Storey ',i+1,end=' ')
#     a=float(input("Force = "))
#     force.append(a)





# height = 10 # Floor height
# storey = 3
# bay = 2
bay_widths = np.array(widths_bay)
#force_list = np.array([10, 8.5, 7, 5, 3.5, 2])
force_list = np.array(force)

bay_widths_cm = np.hstack((0,bay_widths.cumsum()))
heights = height*np.ones(storey)
heights_cm = np.hstack((0,heights.cumsum()))

h_line = storey*bay
v_line = (bay+1)*storey

print('')
print("--------Figure of Frame--------")
def frame():
    fig= plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1,1,1)


    for i in range(storey):
        for j in range(bay):        
            ax.plot([bay_widths_cm[j],bay_widths_cm[j+1]],[heights_cm[i+1], heights_cm[i+1]], color='blue') #Horizontal line

    for i in range(storey):
        for j in range(bay+1):
            ax.plot([bay_widths_cm[j], bay_widths_cm[j]], [heights_cm[i], heights_cm[i+1]], color='blue')
    return fig, ax


fig, ax = frame()
for i in range(storey):
    arl = 5
    ax.annotate(str(force_list[i]), xytext=[-arl, heights_cm[i+1]-0.7],
                 xy=[0, heights_cm[i+1]], 
                 arrowprops=dict(arrowstyle="->", color='red'),
                fontsize=12, size=20)
ax.set_xlim(-bay_widths[0]*0.3, np.max(bay_widths_cm) + bay_widths[0]*0.3)
fig.savefig('Frame_with_lateral_load.png', dpi=300)

shear_unit = force_list[::-1].cumsum()/(bay*2)
csf = np.zeros((storey, bay+1))
csf_shape=np.shape(csf) #shear force for column
cmt = np.zeros((storey, bay+1))#Bending moment for column
bmt = np.zeros((storey, bay))#Bending moment for beam
for i in range(storey):
    csf[i]=shear_unit[i]*2
    csf[i, 0] = shear_unit[i]
    csf[i, -1] = shear_unit[i]

## Shear force diagram for column
print('')
print('---- Shear force diagram for column ----')
fig_f, ax_f = frame()
for i in range(storey):
    for j in range(bay+1):
        ax_f.text(bay_widths_cm[j], heights_cm[i]+(height/2) , csf[::-1][i,j])
        rect = patches.Rectangle((bay_widths_cm[j],heights_cm[i]),
                                 csf[::-1][i,j]*0.4,height,linewidth=1,edgecolor='r',facecolor='green', alpha=0.5)
        ax_f.add_patch(rect)
        #ax_f.text(bay_widths_cm[j]-2, heights_cm[i]+1 , csf[::-1][i,j]*height*0.5)
        #ax_f.text(bay_widths_cm[j]-2, heights_cm[i]+height-1.5 , csf[::-1][i,j]*height*0.5)
ax_f.set_xlim(-bay_widths[0]*0.3, np.max(bay_widths_cm) + bay_widths[0]*0.3)
fig_f.savefig('Shear_force_column.png', dpi=300)


## Bending moment diagram for column
print('')
print('---- Bending moment diagram for column ----')
fig_g, ax_g = frame()
for i in range(storey):
    for j in range(bay+1):
        mnt = csf[::-1][i,j] * height/2
        cmt[i,j]=mnt
        ax_g.text(bay_widths_cm[j], heights_cm[i]+(height/2) , mnt)
        #rect = patches.Rectangle((bay_widths_cm[j],heights_cm[i]),
        #                         csf[::-1][i,j]*0.4,height,linewidth=1,edgecolor='r',facecolor='none')
        ax_g.plot([bay_widths_cm[j], bay_widths_cm[j]+(mnt*0.05),  bay_widths_cm[j]-(mnt*0.05), bay_widths_cm[j]],
                  [heights_cm[i], heights_cm[i], heights_cm[i]+height, heights_cm[i]+height], 
                          color='red')
        #ax_f.add_patch(rect)

ax_g.set_xlim(-bay_widths[0]*0.3, np.max(bay_widths_cm) + bay_widths[0]*0.3)
fig_g.savefig('Moment_diagram_column.png', dpi=300)

## Bending moment diagram for Beam
print('')
print('---- Bending moment diagram for Beam ----')
fig_h, ax_h = frame()
for i in range(storey):
    for j in range(bay):
        bmb = (cmt[i,0]+cmt[i+1,0]) if i+1<=storey-1 else cmt[i,0]
        bmt[i,j]=bmb
        ax_h.text((bay_widths_cm[j]+bay_widths_cm[j+1])/2, heights_cm[i+1] , bmb)
        ax_h.plot([bay_widths_cm[j], bay_widths_cm[j], bay_widths_cm[j+1], bay_widths_cm[j+1]],
                  [heights_cm[i+1], heights_cm[i+1]+(bmb*0.04), heights_cm[i+1]-(bmb*0.04), heights_cm[i+1]], color='red')
        #print ('Storey:',i,'Bay:',j,'Column_shear:',bmb)
        #ax_g.text(bay_widths_cm[j], heights_cm[i]+(height/2) , mnt)


ax_h.set_xlim(-bay_widths[0]*0.3, np.max(bay_widths_cm) + bay_widths[0]*0.3)
fig_h.savefig('Moment_diagram_beam.png', dpi=300)

## Shear force diagram beam

fig_i, ax_i = frame()
for i in range(storey):
    for j in range(bay):
        sf = np.round(bmt[i,j]/(bay_widths[j]*0.5),2) #Shear force
        ax_i.text((bay_widths_cm[j]+bay_widths_cm[j+1])/2, heights_cm[i+1] , sf)
        #ax_i.plot([bay_widths_cm[j], bay_widths_cm[j], bay_widths_cm[j+1], bay_widths_cm[j+1]],
        #          [heights_cm[i+1], heights_cm[i+1]+(sf*0.4), heights_cm[i+1]+(sf*0.4), heights_cm[i+1]], color='red')
        
        
        rect = patches.Rectangle((bay_widths_cm[j],heights_cm[i+1]),
                         bay_widths[j],sf*0.4,linewidth=1,edgecolor='r',facecolor='green', alpha=0.5)
        ax_i.add_patch(rect)
        
print('')
print('---- Shear Force diagram for Beam ----')
ax_i.set_xlim(-bay_widths[0]*0.3, np.max(bay_widths_cm) + bay_widths[0]*0.3)
fig_i.savefig('Shear_diagram_beam.png', dpi=300)