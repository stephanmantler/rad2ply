# --- radio telemetry to Stanford PLY geometry converter ---
#
# written on a stormy day at Grímsfjall (N 64.397 W 17.274)
# stephan mantler / step@stepman.is

import pickle
import numpy as np

# for reprojection
import math

# these only make sense if you run this script from the command line
# input / output file names are ignored when run from within Blender

inputFileName = '20220901_flight2_full_upsample10.pickle'
outputFileName = 'scratch.ply'
# inputFileName = '20220901_flight2_full.pickle'
inputFileName = '20220901_flight2.pickle'

# -- parameters to tune output --
# used directly for command line execution, and set defaults for
# Blender UI

params = {
    # near limit for distance to reflector
    'nearDistanceLimit': 12,
    # far limit for distance to reflector
    'farDistanceLimit': 80,
    # min / max return power to rescale values
    'minReturnPower': -70,
    'maxReturnPower': -40,
    # scaling factor for latitude / longitude and vertical
    'horizontalScale': 0.01,
    'verticalScale': 0.015
}

# local orthographic projection origin
# set to desired value for aligning multiple data sets, or leave np.nan
# to use first data point as origin
# the reprojection code assumes input data is EPSG:4326 (WGS84 lat/lon)
localOriginLat = np.nan
localOriginLon = np.nan

# --- local functions ---


# convert return power to RGBA color
# we just clamp to desired range and convert to greyscale.
def rp_to_color(power, p):
    # -70 to -40 to start with.
    scaled = (power - p['minReturnPower']) / (p['maxReturnPower'] - p['minReturnPower'])  # noqa: E501
    if scaled < 0:
        scaled = 0
    if scaled > 1:
        scaled = 1
    return (scaled * 255, scaled * 255, scaled * 255)


# directly coded orthographic projection.
# Formula and values from EPSG Geomatics Guidance Note Number 7, part 2
def reproject(lat, lon):
    # semi-major axis
    a = 6378137
    # eccentricity of the ellipsoid
    e = 0.081819191
    # projection origin latitude / longitude
    phi0 = localOriginLat * math.pi / 180
    lamda0 = localOriginLon * math.pi / 180
    phi = lat * math.pi / 180
    lamda = lon * math.pi / 180

    e2 = e*e

    cosphi = math.cos(phi)
    sinphi = math.sin(phi)
    cosphi0 = math.cos(phi0)
    sinphi0 = math.sin(phi0)

    # prime verical radius of curvature at phi
    v = a / math.sqrt(1 - e2 * sinphi * sinphi)
    v0 = a / math.sqrt(1 - e2 * sinphi0 * sinphi0)

    # calculate easting and northing
    x = v * cosphi * math.sin(lamda - lamda0)
    y = v * (sinphi * cosphi0 - cosphi * sinphi0 * math.cos(lamda - lamda0))

    y = y + e2 * (v0 * sinphi0 - v * sinphi) * cosphi0

    return ([x, y])


# === main code ===

def readPickle(filename, p):
    global localOriginLat
    global localOriginLon

    print("loading pickle from", filename)

    with open(filename, 'rb') as f:
        data = pickle.load(f)

    # assuming all arrays have the same shape
    # so we pull total sample count from latitude
    # (we should probably check that the pickle has the expected structure)
    timeStepCount = len(data['lat'])

    # figure out distance range
    # (we will disregard power values that are too close / too far away)
    dtr = data['distance_to_reflector']
    distanceCount = len(dtr)

    dtrmin = 0
    dtrmax = distanceCount
    for i in range(distanceCount):
        if dtr[i] < p['nearDistanceLimit']:
            dtrmin = max(i, dtrmin)
        if dtr[i] > p['farDistanceLimit']:
            dtrmax = min(i, dtrmax)

    dtrUsedCount = dtrmax - dtrmin

    print(timeStepCount, "time steps, using range bins", dtrmin, "-", dtrmax)

    plyVertices = []
    plyFaces = []

    print("building geometry ...")

    # makes building geometry and indexing into the vertex array easier
    lastBaseIndex = -1

    for i in range(timeStepCount):

        # early break if we don't have data here
        if np.isnan(data['lat'][i]):
            continue

        if np.isnan(localOriginLat):
            localOriginLat = data['lat'][i]
            localOriginLon = data['lon'][i]

        (x, y) = reproject(data['lat'][i], data['lon'][i])

        # skip samples with invalid positioning
        if np.isnan(x) or np.isnan(y):
            continue

        baseAltitude = data['alt_rel'][i]
        baseIndex = len(plyVertices)
        sxy = p['horizontalScale']
        sz = p['verticalScale']
        for d in range(dtrmin, dtrmax):
            rp = data["return_power"][d, i]
            plyVertices.append({
                "v": (sxy * x, sxy * y, sz * (baseAltitude - dtr[d])),
                "returnPower": rp})

        # only start adding faces if we have at least two columns of vertices
        if lastBaseIndex >= 0:
            for d in range(0, dtrUsedCount-1):
                plyFaces.append([baseIndex + d, baseIndex + d + 1,
                                 lastBaseIndex + d + 1, lastBaseIndex + d])
        lastBaseIndex = baseIndex

    # write everything out
    print(" -", len(plyFaces), "faces generated")
    return (plyVertices, plyFaces)


def writePly(filename, vertices, faces, p):
    print("writing PLY file ...")

    with open(filename, "wb") as file:
        fw = file.write

        fw(b"ply\n")
        fw(b"format ascii 1.0\n")

        fw(b"element vertex %d\n" % len(vertices))
        fw(
            b"property float x\n"
            b"property float y\n"
            b"property float z\n"
        )
        fw(
            b"property uchar red\n"
            b"property uchar green\n"
            b"property uchar blue\n"
            b"property float returnPower\n"
        )

        fw(b"element face %d\n" % len(faces))
        fw(b"property list uchar uint vertex_indices\n")
        fw(b"end_header\n")

        print(" ... vertex data ...")
        for v in vertices:
            fw(b"%.6f %.6f %.6f" % v['v'][:])
            fw(b" %u %u %u" % rp_to_color(v['returnPower'], p)[:])
            fw(b" %.3f" % v['returnPower'])
            fw(b"\n")

        print(" ... face data ...")
        for pf in faces:
            fw(b"%d" % len(pf))
            for index in pf:
                fw(b" %d" % index)
            fw(b"\n")

        print(" ... flushing output buffers ...")
        file.flush()


if __name__ == "__main__":
    (vertices, faces) = readPickle(inputFileName, params)
    writePly(outputFileName, vertices, faces, params)
    print("all done. enjoy!")

# ===== Blender Integration =====

# if this fails we're running outside of Blender
# and there is nothing else to do here.
try:
    import bpy
    from bpy.props import (
        CollectionProperty,
        StringProperty,
        FloatProperty,
    )
    from bpy_extras.io_utils import (
        ImportHelper,
    )
except ImportError:
    exit()

bl_info = {
    "name": "Stanford Glacier Radiometry Import",
    "author": "Stephan Mantler",
    "location": "File > Import/Export",
    "version": (2, 2, 0),
    "blender": (3, 0, 0),
    "description": "Import Stanford Glacier Radiometry data",
    "category": "Import-Export",
}


class ImportStanfordRadiometry(bpy.types.Operator, ImportHelper):
    """Load a PLY geometry file"""
    bl_idname = "import_mesh.stanfordradiometry"  # noqa: F722
    bl_label = "Import Stanford Glacier Radiometry"  # noqa: F722
    bl_options = {'UNDO'}

    files: CollectionProperty(
        name="File Path",  # noqa: F722
        description="File path to the radiometry file",  # noqa: F722
        type=bpy.types.OperatorFileListElement,
    )

    horizontalScale: FloatProperty(
        name="XY Scaling",  # noqa: F722
        description="Horizontal (x/y) multiplier resizes "  # noqa: F722
                    "coordinates to something manageable in Blender",
        default=0.01,
        precision=3,
    )

    verticalScale: FloatProperty(
        name="Z Scaling",  # noqa: F722
        description="Separate vertical (z) multiplier "  # noqa: F722
                    "for vertical exaggeration",
        default=0.015,
        precision=3,
    )

    nearDistanceLimit: FloatProperty(
        name="Near Z Limit",  # noqa: F722
        description="Minimum vertical distance from "  # noqa: F722
                    "receiver to include in geometry",
        default=12,
    )

    farDistanceLimit: FloatProperty(
        name="Far Z Limit",  # noqa: F722
        description="Maximum vertical distance from "  # noqa: F722
                    "receiver to include in geometry",
        default=80,
    )

    directory: StringProperty()

    filename_ext = ".pickle"
    filter_glob: StringProperty(default="*.pickle",  # noqa: F722
                                options={'HIDDEN'})  # noqa: F821

    def execute(self, context):
        import os
        import bmesh

        context.window.cursor_set('WAIT')

        paths = [
            os.path.join(self.directory, name.name)
            for name in self.files
        ]

        if not paths:
            paths.append(self.filepath)

        for ob in bpy.context.selected_objects:
            ob.select_set(False)

        for path in paths:
            name = bpy.path.display_name_from_filepath(path)

            params = {
                # near limit for distance to reflector
                'nearDistanceLimit': self.nearDistanceLimit,
                # far limit for distance to reflector
                'farDistanceLimit': self.farDistanceLimit,
                # scaling factor for latitude / longitude and vertical
                'horizontalScale': self.horizontalScale,
                'verticalScale': self.verticalScale
            }

            print(" -- import:", path)
            (vertices, faces) = readPickle(path, params)
            mesh_data = bpy.data.meshes.new("cube_mesh_data")
            mesh_data.from_pydata([v['v'] for v in vertices], [], faces)
            mesh_data.update()

            bm = bmesh.new()
            bm.from_mesh(mesh_data)
            layer_kind = bm.verts.layers.float.new("returnPower")
            bm.verts.ensure_lookup_table()
            for i, v in enumerate(vertices):
                bm.verts[i][layer_kind] = v['returnPower']

            bm.to_mesh(mesh_data)
            bm.free()

            obj = bpy.data.objects.new(name, mesh_data)
            bpy.context.collection.objects.link(obj)
            bpy.context.view_layer.objects.active = obj
            obj.select_set(True)

            # create suitable material
            material = bpy.data.materials.new(name="mat_"+name)
            material.use_nodes = True

            node_tree = material.node_tree

            # Remove default
            node_tree.nodes.remove(node_tree.nodes.get('Principled BSDF'))
            material_output = node_tree.nodes.get('Material Output')

            attribute_node = node_tree.nodes.new('ShaderNodeAttribute')
            attribute_node.location.x = -600
            attribute_node.attribute_name = "returnPower"

            map_range = node_tree.nodes.new('ShaderNodeMapRange')
            map_range.location.x = -400
            map_range.inputs[1].default_value = -70
            map_range.inputs[2].default_value = -40

            color_ramp = node_tree.nodes.new("ShaderNodeValToRGB")
            cr = color_ramp.color_ramp
            color_ramp.location.x = -200
            cr.elements.new(position=0.14)
            cr.elements.new(position=0.29)
            cr.elements.new(position=0.43)
            cr.elements.new(position=0.57)
            cr.elements.new(position=0.71)
            cr.elements.new(position=0.86)

            cr.elements[0].color = (68/255, 1/255, 84/255, 1)
            cr.elements[1].color = (70/255, 50/255, 127/255, 1)
            cr.elements[2].color = (54/255, 92/255, 141/255, 1)
            cr.elements[3].color = (39/255, 127/255, 142/255, 1)
            cr.elements[4].color = (31/255, 161/255, 135/255, 1)
            cr.elements[5].color = (74/255, 194/255, 109/255, 1)
            cr.elements[6].color = (159/255, 218/255, 58/255, 1)
            cr.elements[7].color = (253/255, 231/255, 3/255, 1)

            diffuse = material.node_tree.nodes.new('ShaderNodeBsdfDiffuse')
            diffuse.location.x = 200

            node_tree.links.new(map_range.inputs[0], attribute_node.outputs[2])
            node_tree.links.new(color_ramp.inputs[0], map_range.outputs[0])
            node_tree.links.new(diffuse.inputs[0], color_ramp.outputs[0])

            # link diffuse shader to material
            node_tree.links.new(material_output.inputs[0], diffuse.outputs[0])

            # set activer material to your new material
            bpy.context.object.active_material = material

            context.window.cursor_set('DEFAULT')

            return {'FINISHED'}


def menu_func_import(self, context):
    self.layout.operator(ImportStanfordRadiometry.bl_idname,
                         text="Stanford Glacier Radiometry (.pickle)")


def register():
    bpy.utils.register_class(ImportStanfordRadiometry)

    bpy.types.TOPBAR_MT_file_import.append(menu_func_import)


def unregister():
    bpy.utils.unregister_class(ImportStanfordRadiometry)

    bpy.types.TOPBAR_MT_file_import.remove(menu_func_import)
