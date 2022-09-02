# --- radio telemetry to Stanford PLY geometry converter ---
#
# written on a stormy day at Grímsfjall (N 64.397 W 17.274)
# stephan mantler / step@stepman.is

import pickle
import numpy as np

inputFileName = '20220901_flight2_full_upsample10.pickle'
outputFileName = 'scratch.ply'
#inputPicke = '20220901_flight2_full.pickle'
#inputPicke = '20220901_flight2.pickle'

# -- parameters to tune output --

# near limit for distance to reflector
distanceLimitNear = 12
# far limit for distance to reflector
distanceLimitFar = 80

# min / max return power to rescale values
minReturnPower = -70
maxReturnPower = -40

# scaling factor for latitude / longitude and vertical
latitudeScale = 10000
longitudeScale = 3000
verticalScale = 0.15

print("loading pickle from",inputFileName)

with open(inputFileName, 'rb') as f:
    data = pickle.load(f)

# assuming all arrays have the same shape so we pull total sample count from latitude
timeStepCount = len(data['lat'])

# figure out distance range
# (we will disregard power values that are too close / too far away)
dtr = data['distance_to_reflector']
distanceCount = len(dtr)

dtrmin = 0
dtrmax = distanceCount
for i in range(distanceCount):
    if dtr[i] < distanceLimitNear:
        dtrmin = max(i, dtrmin)
    if dtr[i] > distanceLimitFar:
        dtrmax = min(i, dtrmax)

dtrUsedCount = dtrmax - dtrmin

print(timeStepCount, "time steps, using range bins", dtrmin, "-", dtrmax)

plyVertices = []
plyFaces = []

lat0 = np.nan
lon0 = np.nan

# convert return power to RGBA color
# we just clamp to desired range and convert to greyscale.
def rp_to_color(power):
    # -70 to -40 to start with.
    scaled = (power - minReturnPower) / (maxReturnPower - minReturnPower)
    if scaled < 0:
        scaled = 0
    if scaled > 1:
        scaled = 1
    return (scaled * 255, scaled * 255, scaled * 255, 255)


print("building geometry ...")

# easy cheat to make building geometry and indexing into the vertex array easier
lastBaseIndex = -1

for i in range(timeStepCount):
    
    # locally rebase all coordinates to the first position we have
    if np.isnan(lat0) and not np.isnan(data['lat'][i]):
        lat0 = data['lat'][i]
        lon0 = data['lon'][i]
        print("local origin set to", lat0, lon0)
    
    lat = (data['lat'][i] - lat0) * latitudeScale
    lon = (data['lon'][i] - lon0) * longitudeScale
    
    # skip samples with invalid positioning
    if np.isnan(lat):
        continue
        
    baseAltitude = data['alt_rel'][i]
    baseIndex = len(plyVertices)

    for d in range(dtrmin, dtrmax):
        rp = data["return_power"][d, i]
        plyVertices.append( { "v": ( lon, lat, verticalScale * (baseAltitude - dtr[d])), "color": rp_to_color(rp) } )

    # only start adding faces if we have at least two columns of vertices
    if lastBaseIndex >= 0:
        for d in range(0, dtrUsedCount-1):            
            plyFaces.append( [ baseIndex + d, baseIndex + d + 1, lastBaseIndex + d + 1, lastBaseIndex + d] )
    lastBaseIndex = baseIndex  

# write everything out
print(" -",len(plyFaces),"faces generated")
print("writing PLY file ...")

with open(outputFileName, "wb") as file:
    fw = file.write

    fw(b"ply\n")
    fw(b"format ascii 1.0\n")
    
    fw(b"element vertex %d\n" % len(plyVertices))
    fw(
        b"property float x\n"
        b"property float y\n"
        b"property float z\n"
    )
    fw(
        b"property uchar red\n"
        b"property uchar green\n"
        b"property uchar blue\n"
        b"property uchar alpha\n"
    )
    
    fw(b"element face %d\n" % len(plyFaces))
    fw(b"property list uchar uint vertex_indices\n")
    fw(b"end_header\n")
    
    # Vertex data
    # ---------------------------
    
    print(" ... vertex data ...")
    for v in plyVertices:
        fw(b"%.6f %.6f %.6f" % v['v'][:])
        fw(b" %u %u %u %u" % v['color'][:])
        fw(b"\n")
    
    # Face data
    # ---------------------------
    
    print(" ... face data ...")
    for pf in plyFaces:
        fw(b"%d" % len(pf))
        for index in pf:
            fw(b" %d" % index)
        fw(b"\n")
        
    print(" ... flushing output buffers ...")
    file.flush()
        
print("all done. enjoy!")