import sys
from importlib.metadata import entry_points

def Viewer(vars, title=None, limits={}, FIPY_VIEWER=None, **kwlimits):
    import os
    errors = []
    attempts = []
    viewers = []

    emptyvars = [var for var in vars if var.mesh.numberOfCells == 0]
    vars = [var for var in vars if var.mesh.numberOfCells > 0]

    if len(emptyvars):
        viewers.append(DummyViewer(vars=emptyvars))

    if FIPY_VIEWER is None:
        enpts = entry_points(group="fipy.viewers")
    else:
        enpts = entry_points(group="fipy.viewers", name=FIPY_VIEWER)

    enpts = sorted(enpts)

    for ep in enpts:
        attempts.append(ep.name)
        try:
            ViewerClass = ep.load()
            while len(vars) > 0:
                viewer = ViewerClass(vars=vars, title=title, limits=limits, **kwlimits)
                for var in viewer.vars:
                    vars.remove(var)
                viewers.append(viewer)
            break
        except Exception as s:
            errors.append("%s: %s" % (ep.name, s))

    return viewers[0]
