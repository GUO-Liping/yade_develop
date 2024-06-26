# encoding: utf-8
"""
Export (not only) geometry to various formats.
"""

from yade.wrapper import *
from yade import utils, Matrix3, Vector3


#textExt===============================================================
def textExt(filename, format='x_y_z_r', comment='', mask=-1, attrs=[]):
	"""Save sphere coordinates and other parameters into a text file in specific format. Non-spherical bodies are silently skipped. Users can add here their own specific format, giving meaningful names. The first file row will contain the format name. Be sure to add the same format specification in ymport.textExt.

	:param string filename: the name of the file, where sphere coordinates will be exported.
	:param string format: the name of output format. Supported 'x_y_z_r'(default), 'x_y_z_r_matId', 'x_y_z_r_attrs' (use proper comment)
	:param string comment: the text, which will be added as a comment at the top of file. If you want to create several lines of text, please use '\\\\n#' for next lines. With 'x_y_z_r_attrs' format, the last (or only) line should consist of column headers of quantities passed as attrs (1 comment word for scalars, 3 comment words for vectors and 9 comment words for matrices)
	:param int mask: export only spheres with the corresponding mask export only spheres with the corresponding mask
	:param [str] attrs: attributes to be exported with 'x_y_z_r_attrs' format. Each str in the list is evaluated for every body exported with body=b (i.e. 'b.state.pos.norm()' would stand for distance of body from coordinate system origin)
	:return: number of spheres which were written.
	:rtype: int
	"""
	O = Omega()

	try:
		out = open(filename, 'w')
	except:
		raise RuntimeError("Problem to write into the file")

	count = 0

	# TODO use output=[] instrad of ''???
	output = ''
	outputVel = ''
	if (format != 'liggghts_in'):
		output = '#format ' + format + '\n'
		if (comment):
			if format == 'x_y_z_r_attrs':
				cmts = comment.split('\n')
				for cmt in cmts[:-1]:
					output += cmt
				output += '# x y z r ' + cmts[-1] + '\n'
			else:
				output += '# ' + comment + '\n'

	minCoord = Vector3.Zero
	maxCoord = Vector3.Zero
	maskNumber = []

	for b in O.bodies:
		try:
			if (isinstance(b.shape, Sphere) and ((mask < 0) or ((mask & b.mask) > 0))):
				if (format == 'x_y_z_r'):
					output += ('%g\t%g\t%g\t%g\n' % (b.state.pos[0], b.state.pos[1], b.state.pos[2], b.shape.radius))
				elif (format == 'x_y_z_r_matId'):
					output += ('%g\t%g\t%g\t%g\t%d\n' % (b.state.pos[0], b.state.pos[1], b.state.pos[2], b.shape.radius, b.material.id))
				elif (format == 'x_y_z_r_attrs'):
					output += ('%g\t%g\t%g\t%g' % (b.state.pos[0], b.state.pos[1], b.state.pos[2], b.shape.radius))
					for cmd in attrs:
						v = eval(cmd)
						if isinstance(v, (int, float)):
							output += '\t%g' % v
						elif isinstance(v, Vector3):
							output += '\t%g\t%g\t%g' % tuple(v[i] for i in range(3))
						elif isinstance(v, Matrix3):
							output += '\t%g' % tuple(v[i] for i in range(9))
					output += '\n'
				elif (format == 'id_x_y_z_r_matId'):
					output += (
					        '%d\t%g\t%g\t%g\t%g\t%d\n' %
					        (b.id, b.state.pos[0], b.state.pos[1], b.state.pos[2], b.shape.radius, b.material.id)
					)
				elif (format == 'jointedPM'):
					output += (
					        '%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n' % (
					                b.id, b.state.onJoint, b.state.joint, b.state.jointNormal1[0], b.state.jointNormal1[1],
					                b.state.jointNormal1[2], b.state.jointNormal2[0], b.state.jointNormal2[1], b.state.jointNormal2[2],
					                b.state.jointNormal3[0], b.state.jointNormal3[1], b.state.jointNormal3[2]
					        )
					)
				elif (format == 'liggghts_in'):
					output += (
					        '%g %g %g %g %g %g %g\n' %
					        (count + 1, b.mask, b.shape.radius, b.material.density, b.state.pos[0], b.state.pos[1], b.state.pos[2])
					)
					outputVel += (
					        '%g %g %g %g %g %g %g\n' % (
					                count + 1, b.state.vel[0], b.state.vel[1], b.state.vel[2], b.state.angVel[0], b.state.angVel[1],
					                b.state.angVel[2]
					        )
					)
				else:
					raise RuntimeError("Please, specify a correct format output!")
				count += 1
				if (count == 1):
					minCoord = b.state.pos - Vector3(b.shape.radius, b.shape.radius, b.shape.radius)
					maxCoord = b.state.pos + Vector3(b.shape.radius, b.shape.radius, b.shape.radius)
				else:
					minCoord = Vector3(
					        min(minCoord[0], b.state.pos[0] - b.shape.radius), min(minCoord[1], b.state.pos[1] - b.shape.radius),
					        min(minCoord[2], b.state.pos[2] - b.shape.radius)
					)
					maxCoord = Vector3(
					        max(maxCoord[0], b.state.pos[0] + b.shape.radius), max(maxCoord[1], b.state.pos[1] + b.shape.radius),
					        max(minCoord[2], b.state.pos[2] + b.shape.radius)
					)
				if b.mask not in maskNumber:
					maskNumber.append(b.mask)
		except AttributeError:
			pass

	if (format == 'liggghts_in'):
		outputHeader = 'LIGGGHTS Description\n\n'
		outputHeader += '%d atoms\n%d atom types\n\n' % (count, len(maskNumber))
		outputHeader += '%g %g xlo xhi\n%g %g ylo yhi\n%g %g zlo zhi\n\n' % (
		        minCoord[0], maxCoord[0], minCoord[1], maxCoord[1], minCoord[2], maxCoord[2]
		)

		output = outputHeader + 'Atoms\n\n' + output + '\nVelocities\n\n' + outputVel
	out.write(output)
	out.close()
	return count


#textClumps===============================================================
def textClumps(filename, format='x_y_z_r_clumpId', comment='', mask=-1):
	"""Save clumps-members into a text file. Non-clumps members are bodies are silently skipped.

	:param string filename: the name of the file, where sphere coordinates will be exported.
	:param string comment: the text, which will be added as a comment at the top of file. If you want to create several lines of text, please use '\\\\n#' for next lines.
	:param int mask: export only spheres with the corresponding mask export only spheres with the corresponding mask
	:return: number of clumps, number of spheres which were written.
	:rtype: int
	"""
	O = Omega()

	try:
		out = open(filename, 'w')
	except:
		raise RuntimeError("Problem to write into the file")

	count = 0
	countClumps = 0
	output = ''
	output = '#format x_y_z_r_clumpId\n'
	if (comment):
		output += '# ' + comment + '\n'

	minCoord = Vector3.Zero
	maxCoord = Vector3.Zero
	maskNumber = []

	for bC in O.bodies:
		if bC.isClump:
			keys = list(bC.shape.members.keys())
			countClumps += 1
			for ii in keys:
				try:
					b = O.bodies[ii]
					if (isinstance(b.shape, Sphere) and ((mask < 0) or ((mask & b.mask) > 0))):
						output += ('%g\t%g\t%g\t%g\t%g\n' % (b.state.pos[0], b.state.pos[1], b.state.pos[2], b.shape.radius, bC.id))
						count += 1
				except AttributeError:
					pass

	out.write(output)
	out.close()
	return countClumps, count


#textPolyhedra===============================================================
def textPolyhedra(fileName, comment='', mask=-1, explanationComment=True, attrs=[]):
	"""Save polyhedra into a text file. Non-polyhedra bodies are silently skipped.

	:param string filename: the name of the output file
	:param string comment: the text, which will be added as a comment at the top of file. If you want to create several lines of text, please use '\\\\n#' for next lines.
	:param int mask: export only polyhedra with the corresponding mask
	:param str explanationComment: inclde explanation of format to the beginning of file
	:return: number of polyhedra which were written.
	:rtype: int
	"""
	count = 0
	f = open(fileName, 'w')
	f.writelines(
	        '# %s\n' % l for l in [
	                'YADE export of polyhedra.',
	                'Each polyhedron export contains first line with id, nuber of vertices and number of surfaces.',
	                'x,y,z coordinates of each vertex follows (each vertex on separate line).',
	                'ids if vertices of individual surfaces follows (numbering from 0, each surface on separate line).',
	                '',
	                'Example of tetrahedron and cube with random ids:',
	                '23 4 4',
	                '0.1 0.2 0.3',
	                '1.3 0.1 -0.1',
	                '-0.2 1.2 0',
	                '0 -0.1 1.5',
	                '0 2 1',
	                '0 3 2',
	                '0 1 3',
	                '1 2 3',
	                '65 8 6',
	                '4 0 0',
	                '5 0 0',
	                '4 1 0',
	                '5 1 0',
	                '4 0 1',
	                '5 0 1',
	                '4 1 1',
	                '5 1 1',
	                '0 2 3 1',
	                '0 1 5 4',
	                '1 3 7 5',
	                '3 2 6 7',
	                '2 0 4 6',
	                '4 5 7 6',
	                '',
	        ]
	)
	if comment:
		f.write('#\n# %s\n' % comment)
	for b in O.bodies:
		if not isinstance(b.shape, Polyhedra) or not mask & b.mask:
			continue
		count += 1
		vertices = [b.state.pos + b.state.ori * v for v in b.shape.v]
		surfaces = b.shape.GetSurfaces()
		strAttrs = ''
		if attrs:
			for cmd in attrs:
				v = eval(cmd)
				if isinstance(v, (int, float)):
					strAttrs += ' %g' % v
				elif isinstance(v, Vector3):
					strAttrs += ' %g %g %g' % tuple(v[i] for i in range(3))
				elif isinstance(v, Matrix3):
					strAttrs += ' %g' % tuple(v[i] for i in range(9))
		f.write('%d %d %d%s\n' % (b.id, len(vertices), len(surfaces), strAttrs))
		f.writelines('%.8e %.8e %.8e\n' % (v[0], v[1], v[2]) for v in vertices)
		f.writelines(' '.join(str(i) for i in surface) + '\n' for surface in surfaces)
	f.close()
	return count


#VTKWriter===============================================================
class VTKWriter(object):
	"""
	USAGE:
	create object vtk_writer = VTKWriter('base_file_name'),
	add to engines PyRunner with command='vtk_writer.snapshot()'
	"""

	def __init__(self, baseName='snapshot', startSnap=0):
		self.snapCount = startSnap
		self.baseName = baseName

	def snapshot(self):
		import xml.dom.minidom
		#import xml.dom.ext # python 2.5 and later

		positions = []
		radii = []

		for b in Omega().bodies:
			if b.mold.name == 'Sphere':
				positions.append(b.phys['se3'][0])
				radii.append(b.mold['radius'])

		# Document and root element
		doc = xml.dom.minidom.Document()
		root_element = doc.createElementNS("VTK", "VTKFile")
		root_element.setAttribute("type", "UnstructuredGrid")
		root_element.setAttribute("version", "0.1")
		root_element.setAttribute("byte_order", "LittleEndian")
		doc.appendChild(root_element)

		# Unstructured grid element
		unstructuredGrid = doc.createElementNS("VTK", "UnstructuredGrid")
		root_element.appendChild(unstructuredGrid)

		# Piece 0 (only one)
		piece = doc.createElementNS("VTK", "Piece")
		piece.setAttribute("NumberOfPoints", str(len(positions)))
		piece.setAttribute("NumberOfCells", "0")
		unstructuredGrid.appendChild(piece)

		### Points ####
		points = doc.createElementNS("VTK", "Points")
		piece.appendChild(points)

		# Point location data
		point_coords = doc.createElementNS("VTK", "DataArray")
		point_coords.setAttribute("type", "Float32")
		point_coords.setAttribute("format", "ascii")
		point_coords.setAttribute("NumberOfComponents", "3")
		points.appendChild(point_coords)

		string = str()
		for x, y, z in positions:
			string += repr(x) + ' ' + repr(y) + ' ' + repr(z) + ' '
		point_coords_data = doc.createTextNode(string)
		point_coords.appendChild(point_coords_data)

		#### Cells ####
		cells = doc.createElementNS("VTK", "Cells")
		piece.appendChild(cells)

		# Cell locations
		cell_connectivity = doc.createElementNS("VTK", "DataArray")
		cell_connectivity.setAttribute("type", "Int32")
		cell_connectivity.setAttribute("Name", "connectivity")
		cell_connectivity.setAttribute("format", "ascii")
		cells.appendChild(cell_connectivity)

		# Cell location data
		connectivity = doc.createTextNode("0")
		cell_connectivity.appendChild(connectivity)

		cell_offsets = doc.createElementNS("VTK", "DataArray")
		cell_offsets.setAttribute("type", "Int32")
		cell_offsets.setAttribute("Name", "offsets")
		cell_offsets.setAttribute("format", "ascii")
		cells.appendChild(cell_offsets)
		offsets = doc.createTextNode("0")
		cell_offsets.appendChild(offsets)

		cell_types = doc.createElementNS("VTK", "DataArray")
		cell_types.setAttribute("type", "UInt8")
		cell_types.setAttribute("Name", "types")
		cell_types.setAttribute("format", "ascii")
		cells.appendChild(cell_types)
		types = doc.createTextNode("1")
		cell_types.appendChild(types)

		#### Data at Points ####
		point_data = doc.createElementNS("VTK", "PointData")
		piece.appendChild(point_data)

		# Particle radii
		if len(radii) > 0:
			radiiNode = doc.createElementNS("VTK", "DataArray")
			radiiNode.setAttribute("Name", "radii")
			radiiNode.setAttribute("type", "Float32")
			radiiNode.setAttribute("format", "ascii")
			point_data.appendChild(radiiNode)

			string = str()
			for r in radii:
				string += repr(r) + ' '
			radiiData = doc.createTextNode(string)
			radiiNode.appendChild(radiiData)

		#### Cell data (dummy) ####
		cell_data = doc.createElementNS("VTK", "CellData")
		piece.appendChild(cell_data)

		# Write to file and exit
		outFile = open(self.baseName + '%08d' % self.snapCount + '.vtu', 'w')
		#		xml.dom.ext.PrettyPrint(doc, file)
		doc.writexml(outFile, newl='\n')
		outFile.close()
		self.snapCount += 1


#text===============================================================
def text(filename, mask=-1):
	"""Save sphere coordinates into a text file; the format of the line is: x y z r. Non-spherical bodies are silently skipped. Example added to examples/regular-sphere-pack/regular-sphere-pack.py

	:param string filename: the name of the file, where sphere coordinates will be exported.
	:param int mask: export only spheres with the corresponding mask
	:return: number of spheres which were written.
	:rtype: int
	"""
	return (textExt(filename=filename, format='x_y_z_r', mask=mask))


#VTKExporter===============================================================


class VTKExporter(object):
	"""Class for exporting data to `VTK Simple Legacy File <https://www.vtk.org/VTK/img/file-formats.pdf>`_ (for example if, for some reason, you are not able to use :yref:`VTKRecorder`).
	Supported export of:
	
	* spheres
	* facets
	* polyhedra
	* PotentialBlocks
	* interactions
	* contact points
	* periodic cell
	
	Usage:

	* create object ``vtkExporter = VTKExporter('baseFileName')``,
	* add to ``O.engines`` a ``PyRunner`` with ``command='vtkExporter.exportSomething(...)'``
	* alternatively, just use ``vtkExporter.exportSomething(...)`` at the end of the script for instance
	
	Example: :ysrc:`examples/test/vtk-exporter/vtkExporter.py`, :ysrc:`examples/test/unv-read/unvReadVTKExport.py`.
	
	:param string baseName: name of the exported files. The files would be named, e.g., ``baseName-spheres-snapNb.vtk`` or ``baseName-facets-snapNb.vtk``
	:param int startSnap: the numbering of files will start form ``startSnap``
	"""

	# TODO comments
	def __init__(self, baseName, startSnap=0):
		self.spheresSnapCount = startSnap
		self.facetsSnapCount = startSnap
		self.intrsSnapCount = startSnap
		self.polyhedraSnapCount = startSnap
		self.PotentialBlocksSnapCount = startSnap
		self.contactPointsSnapCount = startSnap
		self.baseName = baseName

	# auxiliary functions
	def _checkWhatArgumentIsDict(self, what, funName, whatName="what"):
		"""An auxiliary function, to be deleted when the 'deprecation period' is over"""
		if isinstance(what, (tuple, list)):
			raise DeprecationWarning("{}: '{}' argument is no longer list/tuple, but dict".format(funName, whatName))
		assert isinstance(what, dict)

	def _warn(self, msg):
		print("Warning (yade.export.VTKExporter): " + msg)

	def _error(self, msg):
		print("ERROR (yade.export.VTKExporter): " + msg)

	def _getBodies(self, ids, type):
		allIds = False
		if isinstance(ids, str) and ids.lower() == 'all':
			ids = range(len(O.bodies))
			allIds = True
		bodies = []
		for i in ids:
			b = O.bodies[i]
			if not b:
				continue
			if not isinstance(b.shape, type):
				if not allIds:
					self._warn("body %d is not of type %s" % (i, type))
				continue
			bodies.append(b)
			if not bodies:
				self._warn("no bodies...")
		return bodies

	def _getInteractions(self, ids):
		if isinstance(ids, str) and ids.lower() == 'all':
			ids = [(i.id1, i.id2) for i in O.interactions]
		intrs = [(i, j) for i, j in ids]
		if not intrs:
			self._warn("no interactions ...")
		return intrs

	def exportSpheres(self, ids='all', what={}, comment="comment", numLabel=None, useRef=False):
		"""exports spheres (positions and radius) and defined properties.
		
		:param [int]|"all" ids: if "all", then export all spheres, otherwise only spheres from integer list
		:param dictionary what: which additional quantities (other than the position and the radius) to export. parameter is name->command dictionary. Name is string under which it is save to vtk, command is string to evaluate. Note that the bodies are labeled as b in this function. Scalar, vector and tensor variables are supported. For example, to export velocity (with name particleVelocity) and the distance form point (0,0,0) (named as dist) you should write: ``what=dict(particleVelocity='b.state.vel',dist='b.state.pos.norm()', ... )``
		:param string comment: comment to add to vtk file
		:param int numLabel: number of file (e.g. time step), if unspecified, the last used value + 1 will be used
		:param bool useRef: if False (default), use current position of the spheres for export, use reference position otherwise
		"""
		# get list of bodies to export
		bodies = self._getBodies(ids, Sphere)
		if not bodies:
			return
		nBodies = len(bodies)
		# output file
		fName = self.baseName + '-spheres-%08d' % (numLabel if numLabel else self.spheresSnapCount) + '.vtk'
		outFile = open(fName, 'w')
		# head
		outFile.write("# vtk DataFile Version 3.0.\n%s\nASCII\n\nDATASET POLYDATA\nPOINTS %d double\n" % (comment, nBodies))
		# write position of spheres
		for b in bodies:
			pos = b.state.refPos if useRef else b.state.pos if not O.periodic else O.cell.wrap(b.state.pos)
			outFile.write("%g %g %g\n" % (pos[0], pos[1], pos[2]))
		# write radius
		outFile.write("\nPOINT_DATA %d\nSCALARS radius double 1\nLOOKUP_TABLE default\n" % (nBodies))
		for b in bodies:
			outFile.write("%g\n" % (b.shape.radius))
		# checks what argument
		self._checkWhatArgumentIsDict(what, "exportSpheres")
		# write additional data from 'what' param
		for name, command in what.items():  # for each name...
			test = eval(command)  # ... eval one example to see what type (float, Vector3, Matrix3) the result is ...
			# ... and write appropriate header line and loop over all bodies and write appropriate vtk line(s)
			if isinstance(test, Matrix3):
				outFile.write("\nTENSORS %s double\n" % (name))
				for b in bodies:
					t = eval(command)
					outFile.write(
					        "%g %g %g\n%g %g %g\n%g %g %g\n\n" %
					        (t[0, 0], t[0, 1], t[0, 2], t[1, 0], t[1, 1], t[1, 2], t[2, 0], t[2, 1], t[2, 2])
					)
			elif isinstance(test, Vector3):
				outFile.write("\nVECTORS %s double\n" % (name))
				for b in bodies:
					v = eval(command)
					outFile.write("%g %g %g\n" % (v[0], v[1], v[2]))
			elif isinstance(test, (int, float)):
				outFile.write("\nSCALARS %s double 1\nLOOKUP_TABLE default\n" % (name))
				for b in bodies:
					outFile.write("%g\n" % (eval(command)))
			else:
				self._warn("exportSpheres: wrong 'what' parameter, vtk output might be corrupted'")
		outFile.close()
		self.spheresSnapCount += 1

	def exportFacets(self, ids='all', what={}, comment="comment", numLabel=None):
		"""
		exports facets (positions) and defined properties. Facets are exported with multiplicated nodes
		
		:param [int]|"all" ids: if "all", then export all facets, otherwise only facets from integer list
		:param dictionary what: see :meth:`exportSpheres`
		:param string comment: comment to add to vtk file
		:param int numLabel: number of file (e.g. time step), if unspecified, the last used value + 1 will be used
		"""
		# get list of bodies to export
		bodies = self._getBodies(ids, Facet)
		if not bodies:
			return
		nBodies = len(bodies)
		# output file
		fName = self.baseName + '-facets-%08d' % (numLabel if numLabel else self.facetsSnapCount) + '.vtk'
		outFile = open(fName, 'w')
		# head
		outFile.write("# vtk DataFile Version 3.0.\n%s\nASCII\n\nDATASET POLYDATA\nPOINTS %d double\n" % (comment, 3 * nBodies))
		# write vertices
		for b in bodies:
			p = b.state.pos
			o = b.state.ori
			s = b.shape
			pt1 = p + o * s.vertices[0]
			pt2 = p + o * s.vertices[1]
			pt3 = p + o * s.vertices[2]
			outFile.write("%g %g %g\n" % (pt1[0], pt1[1], pt1[2]))
			outFile.write("%g %g %g\n" % (pt2[0], pt2[1], pt2[2]))
			outFile.write("%g %g %g\n" % (pt3[0], pt3[1], pt3[2]))
		# write facets
		outFile.write("\nPOLYGONS %d %d\n" % (nBodies, 4 * nBodies))
		i = 0
		for b in bodies:
			outFile.write("3 %d %d %d\n" % (i, i + 1, i + 2))
			i += 3
		# checks what argument
		self._checkWhatArgumentIsDict(what, "exportFacets")
		# write additional data from 'what' param
		if what:
			outFile.write("\nCELL_DATA %d" % (nBodies))
		# see exportSpheres for explanation of this code block
		for name, command in what.items():
			test = eval(command)
			if isinstance(test, Matrix3):
				outFile.write("\nTENSORS %s double\n" % (name))
				for b in bodies:
					t = eval(command)
					outFile.write(
					        "%g %g %g\n%g %g %g\n%g %g %g\n\n" %
					        (t[0, 0], t[0, 1], t[0, 2], t[1, 0], t[1, 1], t[1, 2], t[2, 0], t[2, 1], t[2, 2])
					)
			if isinstance(test, Vector3):
				outFile.write("\nVECTORS %s double\n" % (name))
				for b in bodies:
					v = eval(command)
					outFile.write("%g %g %g\n" % (v[0], v[1], v[2]))
			else:
				outFile.write("\nSCALARS %s double 1\nLOOKUP_TABLE default\n" % (name))
				for b in bodies:
					outFile.write("%g\n" % (eval(command)))
		outFile.close()
		self.facetsSnapCount += 1

	def exportFacetsAsMesh(self, ids='all', connectivityTable=None, what={}, comment="comment", numLabel=None):
		"""
		exports facets (positions) and defined properties. Facets are exported as mesh (not with multiplicated nodes). Therefore additional parameters connectivityTable is needed
		
		:param [int]|"all" ids: if "all", then export all facets, otherwise only facets from integer list
		:param dictionary what: see :meth:`exportSpheres`
		:param string comment: comment to add to vtk file
		:param int numLabel: number of file (e.g. time step), if unspecified, the last used value + 1 will be used
		:param [(float,float,float)|Vector3] nodes: list of coordinates of nodes
		:param [(int,int,int)] connectivityTable: list of node ids of individual elements (facets)
		"""
		# get list of bodies to export
		bodies = self._getBodies(ids, Facet)
		ids = [b.id for b in bodies]
		if not bodies:
			return
		nBodies = len(bodies)
		if connectivityTable is None:
			self._error("'connectivityTable' not specified")
			return
		if nBodies != len(connectivityTable):
			self._error("length of 'connectivityTable' does not match length of 'ids', no export")
			return
		# nodes
		nodes = [Vector3.Zero for i in range(max(max(e) for e in connectivityTable) + 1)]
		for id, e in zip(ids, connectivityTable):
			b = bodies[id]
			p = b.state.pos
			o = b.state.ori
			s = b.shape
			pt1 = p + o * s.vertices[0]
			pt2 = p + o * s.vertices[1]
			pt3 = p + o * s.vertices[2]
			nodes[e[0]] = pt1
			nodes[e[1]] = pt2
			nodes[e[2]] = pt3
		# output file
		fName = self.baseName + '-facets-%08d' % (numLabel if numLabel else self.facetsSnapCount) + '.vtk'
		outFile = open(fName, 'w')
		# head
		outFile.write("# vtk DataFile Version 3.0.\n%s\nASCII\n\nDATASET POLYDATA\nPOINTS %d double\n" % (comment, len(nodes)))
		# write vertices
		for node in nodes:
			outFile.write("%g %g %g\n" % (node[0], node[1], node[2]))
		# write facets
		outFile.write("\nPOLYGONS %d %d\n" % (len(connectivityTable), 4 * len(connectivityTable)))
		for e in connectivityTable:
			outFile.write("3 %d %d %d\n" % e)
		# checks what argument
		self._checkWhatArgumentIsDict(what, "exportFacetsAsMesh")
		# write additional data from 'what' param
		if what:
			outFile.write("\nCELL_DATA %d" % (nBodies))
		# see exportSpheres for explanation of this code block
		for name, command in what.items():
			test = eval(command)
			if isinstance(test, Matrix3):
				outFile.write("\nTENSORS %s double\n" % (name))
				for b in bodies:
					t = eval(command)
					outFile.write(
					        "%g %g %g\n%g %g %g\n%g %g %g\n\n" %
					        (t[0, 0], t[0, 1], t[0, 2], t[1, 0], t[1, 1], t[1, 2], t[2, 0], t[2, 1], t[2, 2])
					)
			if isinstance(test, Vector3):
				outFile.write("\nVECTORS %s double\n" % (name))
				for b in bodies:
					v = eval(command)
					outFile.write("%g %g %g\n" % (v[0], v[1], v[2]))
			else:
				outFile.write("\nSCALARS %s double 1\nLOOKUP_TABLE default\n" % (name))
				for b in bodies:
					outFile.write("%g\n" % (eval(command)))
		outFile.close()
		self.facetsSnapCount += 1

	def exportInteractions(self, ids='all', what={}, verticesWhat={}, comment="comment", numLabel=None, useRef=False):
		"""exports interactions and defined properties.
		
		:param [(int,int)]|"all" ids: if "all", then export all interactions, otherwise only interactions from (int,int) list
		:param dictionary what: what to export. parameter is a name->command dictionary. Name is string under which it is saved to vtk, command is string to evaluate. Note that the interactions are labeled as i in this function. Scalar, vector and tensor variables are supported. For example, to export the stiffness difference (named as ``dStiff``) from a certain value (1e9) you should write: ``what=dict(dStiff='i.phys.kn-1e9', ... )``
		:param dictionary verticesWhat: what to export on connected bodies. Bodies are labeled as ``b`` (or ``b1`` and ``b2`` if you need to treat both bodies differently)
		:param string comment: comment to add to vtk file
		:param int numLabel: number of file (e.g. time step), if unspecified, the last used value + 1 will be used
		:param bool useRef: if False (default), use current position of the bodies for export, use reference position otherwise
		"""
		# get list of interactions to export
		intrs = self._getInteractions(ids)
		if not intrs:
			return
		nIntrs = len(intrs)
		# output file
		fName = self.baseName + '-intrs-%08d' % (numLabel if numLabel else self.intrsSnapCount) + '.vtk'
		outFile = open(fName, 'w')
		# head
		outFile.write("# vtk DataFile Version 3.0.\n%s\nASCII\n\nDATASET POLYDATA\nPOINTS %d double\n" % (comment, 2 * nIntrs))
		# write coords of intrs bodies (also taking into account possible periodicity
		for ii, jj in intrs:
			i = O.interactions[ii, jj]
			pos = O.bodies[ii].state.refPos if useRef else O.bodies[ii].state.pos
			outFile.write("%g %g %g\n" % (pos[0], pos[1], pos[2]))
			pos = (O.bodies[jj].state.refPos if useRef else O.bodies[jj].state.pos) + (O.cell.hSize * i.cellDist if O.periodic else Vector3.Zero)
			outFile.write("%g %g %g\n" % (pos[0], pos[1], pos[2]))
		# write interactions as lines
		outFile.write("LINES %d %d\n" % (nIntrs, 3 * nIntrs))
		for j, i in enumerate(intrs):
			outFile.write("2 %d %d\n" % (2 * j, 2 * j + 1))
		# checks what argument
		self._checkWhatArgumentIsDict(what, "exportInteractions")
		# write additional data from 'what' param
		if what:
			outFile.write("\nCELL_DATA %d\n" % (nIntrs))
		for i in O.interactions:
			if i.isReal:
				break
		# see exportSpheres for explanation of this code block
		for name, command in what.items():
			test = eval(command)
			if isinstance(test, Matrix3):
				outFile.write("\nTENSORS %s double\n" % (name))
				for ii, jj in intrs:
					i = O.interactions[ii, jj]
					t = eval(command)
					outFile.write(
					        "%g %g %g\n%g %g %g\n%g %g %g\n\n" %
					        (t[0, 0], t[0, 1], t[0, 2], t[1, 0], t[1, 1], t[1, 2], t[2, 0], t[2, 1], t[2, 2])
					)
			elif isinstance(test, Vector3):
				outFile.write("\nVECTORS %s double\n" % (name))
				for ii, jj in intrs:
					i = O.interactions[ii, jj]
					v = eval(command)
					outFile.write("%g %g %g\n" % (v[0], v[1], v[2]))
			elif isinstance(test, (int, float)):
				outFile.write("\nSCALARS %s double 1\nLOOKUP_TABLE default\n" % (name))
				for ii, jj in intrs:
					i = O.interactions[ii, jj]
					outFile.write("%g\n" % (eval(command)))
			else:
				self._warn("exportInteractions: wrong 'what' parameter, vtk output might be corrupted")
		# checks verticesWhat argument
		self._checkWhatArgumentIsDict(verticesWhat, "exportInteractions", "verticesWhat")
		# write additional data of bodies
		if verticesWhat:
			outFile.write("\nPOINT_DATA %d\n" % (2 * nIntrs))
			b = b1 = b2 = O.bodies[0]
		# see exportSpheres for explanation of this code block
		for name, vWhat in verticesWhat.items():
			lw = len(vWhat)
			if lw == 1:
				command = vWhat
				test = eval(command)
			elif lw == 2:
				command1, command2 = vWhat
				test = eval(command1)
			if isinstance(test, Matrix3):
				outFile.write("\nTENSORS %s double\n" % (name))
				for ii, jj in intrs:
					i = O.interactions[ii, jj]
					b1 = O.bodies[ii]
					b2 = O.bodies[jj]
					if lw == 2:
						for b in (b1, b2):
							t = eval(command)
							outFile.write(
							        "%g %g %g\n%g %g %g\n%g %g %g\n\n" %
							        (t[0, 0], t[0, 1], t[0, 2], t[1, 0], t[1, 1], t[1, 2], t[2, 0], t[2, 1], t[2, 2])
							)
					elif lw == 3:
						t1 = eval(command1)
						t2 = eval(command2)
						outFile.write(
						        "%g %g %g\n%g %g %g\n%g %g %g\n\n" %
						        (t1[0, 0], t1[0, 1], t1[0, 2], t1[1, 0], t1[1, 1], t1[1, 2], t1[2, 0], t1[2, 1], t1[2, 2])
						)
					outFile.write(
					        "%g %g %g\n%g %g %g\n%g %g %g\n\n" %
					        (t2[0, 0], t2[0, 1], t2[0, 2], t2[1, 0], t2[1, 1], t2[1, 2], t2[2, 0], t2[2, 1], t2[2, 2])
					)
			elif isinstance(test, Vector3):
				outFile.write("\nVECTORS %s double\n" % (name))
				for ii, jj in intrs:
					i = O.interactions[ii, jj]
					b1 = O.bodies[ii]
					b2 = O.bodies[jj]
					if lw == 2:
						for b in (b1, b2):
							v = eval(command)
							outFile.write("%g %g %g\n" % (v[0], v[1], v[2]))
					elif lw == 3:
						v1 = eval(command1)
						v2 = eval(command2)
						outFile.write("%g %g %g\n" % (v1[0], v1[1], v1[2]))
						outFile.write("%g %g %g\n" % (v2[0], v2[1], v2[2]))
			elif isinstance(test, (int, float)):
				outFile.write("\nSCALARS %s double 1\nLOOKUP_TABLE default\n" % (name))
				for ii, jj in intrs:
					i = O.interactions[ii, jj]
					b1 = O.bodies[ii]
					b2 = O.bodies[jj]
					if lw == 2:
						for b in (b1, b2):
							outFile.write("%g\n" % (eval(command)))
					elif lw == 3:
						outFile.write("%g\n" % (eval(command1)))
						outFile.write("%g\n" % (eval(command2)))
			else:
				self._warn("exportInteractions: wrong 'what' parameter, vtk output might be corrupted")
		outFile.close()
		self.intrsSnapCount += 1

	def exportContactPoints(self, ids='all', what={}, useRef={}, comment="comment", numLabel=None):
		"""exports contact points (CPs) and defined properties.
		
		:param [(int,int)] ids: see :meth:`exportInteractions`
		:param dictionary what: see :meth:`exportInteractions`
		:param bool useRef: see :meth:`exportInteractions`
		:param string comment: comment to add to vtk file
		:param int numLabel: number of file (e.g. time step), if unspecified, the last used value + 1 will be used
		"""
		# get list of interactions to export
		if useRef:
			useRef = dict(((i.id1, i.id2), v) for i, v in useRef.items())
			intrs = list(useRef.keys())
		else:
			intrs = self._getInteractions(ids)
		if not intrs:
			return
		nIntrs = len(intrs)
		# output file
		fName = self.baseName + '-cps-%08d' % (numLabel if numLabel else self.contactPointsSnapCount) + '.vtk'
		outFile = open(fName, 'w')
		# head
		outFile.write("# vtk DataFile Version 3.0.\n%s\nASCII\n\nDATASET POLYDATA\nPOINTS %d double\n" % (comment, nIntrs))
		# write coords of contact points
		for ii, jj in intrs:
			if useRef:
				pos = useRef[(ii, jj)]
			else:
				i = O.interactions[ii, jj]
				pos = i.geom.contactPoint
			outFile.write("%g %g %g\n" % (pos[0], pos[1], pos[2]))
		# checks what argument
		self._checkWhatArgumentIsDict(what, "exportContactPoints")
		# see exportSpheres for explanation of this code block
		if what:
			outFile.write("\nPOINT_DATA %d\n" % (nIntrs))
			for i in O.interactions:
				if i.isReal:
					break
		for name, command in what.items():
			test = eval(command)
			if isinstance(test, Matrix3):
				outFile.write("\nTENSORS %s double\n" % (name))
				for ii, jj in intrs:
					try:
						i = O.interactions[ii, jj]
						t = eval(command)
					except IndexError:
						t = Matrix3.Zero  # TODO?
					outFile.write(
					        "%g %g %g\n%g %g %g\n%g %g %g\n\n" %
					        (t[0, 0], t[0, 1], t[0, 2], t[1, 0], t[1, 1], t[1, 2], t[2, 0], t[2, 1], t[2, 2])
					)
			elif isinstance(test, Vector3):
				outFile.write("\nVECTORS %s double\n" % (name))
				for ii, jj in intrs:
					try:
						i = O.interactions[ii, jj]
						v = eval(command)
					except IndexError:
						v = Vector3.Zero  # TODO?
					outFile.write("%g %g %g\n" % (v[0], v[1], v[2]))
			elif isinstance(test, (int, float)):
				outFile.write("\nSCALARS %s double 1\nLOOKUP_TABLE default\n" % (name))
				for ii, jj in intrs:
					try:
						i = O.interactions[ii, jj]
						f = eval(command)
					except IndexError:
						f = 0.  # TODO?
					outFile.write("%g\n" % (f))
			else:
				self._warn("exportContacPoints: wrong 'what' parameter, vtk output might be corrupted'")
		outFile.close()
		self.contactPointsSnapCount += 1

	def exportPeriodicCell(self, comment="comment", numLabel=None):
		"""exports the :yref:`Cell` geometry for periodic simulations.
		
		:param string comment: comment to add to vtk file
		:param int numLabel: number of file (e.g. time step), if unspecified, the last used value + 1 will be used
		"""
		if not O.periodic:
			self._warn("exportPeriodicCell: scene is not periodic, no export...")
			return
		hSize = O.cell.hSize
		fName = self.baseName + '-periCell-%08d' % (numLabel if numLabel else self.intrsSnapCount) + '.vtk'
		outFile = open(fName, 'w')
		outFile.write("# vtk DataFile Version 3.0.\n%s\nASCII\n\nDATASET UNSTRUCTURED_GRID\nPOINTS 8 double\n" % (comment))
		vertices = [
		        hSize * Vector3(0, 0, 1),
		        hSize * Vector3(0, 1, 1),
		        hSize * Vector3(1, 1, 1),
		        hSize * Vector3(1, 0, 1),
		        hSize * Vector3(0, 0, 0),
		        hSize * Vector3(0, 1, 0),
		        hSize * Vector3(1, 1, 0),
		        hSize * Vector3(1, 0, 0),
		]
		for v in vertices:
			outFile.write('%g %g %g\n' % (v[0], v[1], v[2]))
		outFile.write('\nCELLS 1 9\n')
		outFile.write('8 0 1 2 3 4 5 6 7\n')
		outFile.write('\nCELL_TYPES 1\n12\n')
		outFile.close()

	def exportPolyhedra(self, ids='all', what={}, comment="comment", numLabel=None, useRef=False):
		"""Exports polyhedrons and defined properties.
		
		:param ids: if "all", then export all polyhedrons, otherwise only polyhedrons from integer list
		:type ids: [int] | "all"
		:param dictionary what: which additional quantities (in addition to the positions) to export. parameter is name->command dictionary. Name is string under which it is saved to vtk, command is string to evaluate. Note that the bodies are labeled as b in this function. Scalar, vector and tensor variables are supported. For example, to export velocity (named as particleVelocity) and the distance from point (0,0,0) (named as dist) you should write: ``what=dict(particleVelocity='b.state.vel',dist='b.state.pos.norm()', ... )``
		:param string comment: comment to add to vtk file
		:param int numLabel: number of file (e.g. time step), if unspecified, the last used value + 1 will be used
		"""
		# TODO useRef?
		# get list of bodies to export
		bodies = self._getBodies(ids, Polyhedra)  # TODO
		if not bodies:
			return
		# number of vertices
		nVertices = sum(len(b.shape.v) for b in bodies)
		# export polyherda as a set of triangle faces
		bodyFaces = []
		for b in bodies:
			ff = []
			f = b.shape.GetSurfaceTriangulation()
			for i in range(len(f) // 3):
				ff.append([f[3 * i + j] for j in (0, 1, 2)])
			bodyFaces.append(ff)
		# output file
		nFaces = sum(len(f) for f in bodyFaces)
		fName = self.baseName + '-polyhedra-%08d' % (numLabel if numLabel else self.polyhedraSnapCount) + '.vtk'
		outFile = open(fName, 'w')
		# head
		outFile.write("# vtk DataFile Version 3.0.\n%s\nASCII\n\nDATASET POLYDATA\nPOINTS %d double\n" % (comment, nVertices))
		# write position of vertices
		if useRef:
			dspls = []
			for b in bodies:
				bPos = b.state.pos
				bOri = b.state.ori
				brPos = b.state.refPos
				brOri = b.state.refOri
				for v in b.shape.v:
					rPos = brPos + brOri * v
					pos = bPos + bOri * v
					outFile.write("%g %g %g\n" % (rPos[0], rPos[1], rPos[2]))
					dspls.append(pos - rPos)
		else:
			for b in bodies:
				bPos = b.state.pos
				bOri = b.state.ori
				for v in b.shape.v:
					pos = bPos + bOri * v
					outFile.write("%g %g %g\n" % (pos[0], pos[1], pos[2]))
		# write triangle faces
		outFile.write("\nPOLYGONS %d %d\n" % (nFaces, 4 * nFaces))
		j = 0
		for i, b in enumerate(bodies):
			faces = bodyFaces[i]
			for face in faces:
				t = tuple([j + ii for ii in face])
				outFile.write("3 %d %d %d\n" % t)
			j += len(b.shape.v)
		# checks what argument
		self._checkWhatArgumentIsDict(what, "exportPolyhedra")
		# write additional data from 'what' param
		if useRef:
			outFile.write("\nPOINT_DATA %d\n" % (len(dspls)))
			outFile.write("\nVECTORS displacement double\n")
			for v in dspls:
				outFile.write("%g %g %g\n" % (v[0], v[1], v[2]))
		if what:
			outFile.write("\nCELL_DATA %d" % (nFaces))
		# see exportSpheres for explanation of this code block
		for name, command in what.items():
			test = eval(command)
			if isinstance(test, Matrix3):
				outFile.write("\nTENSORS %s double\n" % (name))
				for i, b in enumerate(bodies):
					t = eval(command)
					for f in bodyFaces[i]:
						outFile.write(
						        "%g %g %g\n%g %g %g\n%g %g %g\n\n" %
						        (t[0, 0], t[0, 1], t[0, 2], t[1, 0], t[1, 1], t[1, 2], t[2, 0], t[2, 1], t[2, 2])
						)
			elif isinstance(test, Vector3):
				outFile.write("\nVECTORS %s double\n" % (name))
				for i, b in enumerate(bodies):
					v = eval(command)
					for f in bodyFaces[i]:
						outFile.write("%g %g %g\n" % (v[0], v[1], v[2]))
			elif isinstance(test, (int, float)):
				outFile.write("\nSCALARS %s double 1\nLOOKUP_TABLE default\n" % (name))
				for i, b in enumerate(bodies):
					e = eval(command)
					for f in bodyFaces[i]:
						outFile.write("%g\n" % e)
			else:
				self._warn("exportPolyhedra: wrong 'what' parameter, vtk output might be corrupted")
		outFile.close()
		self.polyhedraSnapCount += 1

	def exportPotentialBlocks(self, ids='all', what={}, comment="comment", numLabel=None, useRef=False):
		"""Exports Potential Blocks and defined properties.
		
		:param ids: if "all", then export all Potential Blocks, otherwise only Potential Blocks from integer list
		:type ids: [int] | "all"
		:param dictionary what: which additional quantities (in addition to the positions) to export. parameter is name->command dictionary. Name is string under which it is saved to vtk, command is string to evaluate. Note that the bodies are labeled as b in this function. Scalar, vector and tensor variables are supported. For example, to export velocity (named as particleVelocity) and the distance from point (0,0,0) (named as dist) you should write: ``what=dict(particleVelocity='b.state.vel',dist='b.state.pos.norm()', ... )``
		:param string comment: comment to add to vtk file
		:param int numLabel: number of file (e.g. time step), if unspecified, the last used value + 1 will be used
		"""
		# TODO useRef?
		# get list of bodies to export
		bodies = self._getBodies(ids, PotentialBlock)  # TODO
		if not bodies:
			return
		# number of vertices
		nVertices = sum(len(b.shape.vertices) for b in bodies)
		# export triangulation of particle face
		bodyFaces = []
		for b in bodies:
			ff = []
			for f in b.shape.connectivity:
				for i in range(1, len(f) - 1):
					ff.append([f[0], f[i], f[i + 1]])
			bodyFaces.append(ff)
		# output file
		nFaces = sum(len(f) for f in bodyFaces)
		fName = self.baseName + '-PotentialBlocks-%08d' % (numLabel if numLabel else self.PotentialBlocksSnapCount) + '.vtk'
		outFile = open(fName, 'w')
		# head
		outFile.write("# vtk DataFile Version 3.0.\n%s\nASCII\n\nDATASET POLYDATA\nPOINTS %d double\n" % (comment, nVertices))
		# write position of vertices
		if useRef:
			dspls = []
			for b in bodies:
				bPos = b.state.pos
				bOri = b.state.ori
				brPos = b.state.refPos
				brOri = b.state.refOri
				for v in b.shape.vertices:
					rPos = brPos + brOri * v
					pos = bPos + bOri * v
					outFile.write("%g %g %g\n" % (rPos[0], rPos[1], rPos[2]))
					dspls.append(pos - rPos)
		else:
			for b in bodies:
				bPos = b.state.pos
				bOri = b.state.ori
				for v in b.shape.vertices:
					pos = bPos + bOri * v
					outFile.write("%g %g %g\n" % (pos[0], pos[1], pos[2]))
		# write triangle faces
		outFile.write("\nPOLYGONS %d %d\n" % (nFaces, 4 * nFaces))
		j = 0
		for i, b in enumerate(bodies):
			faces = bodyFaces[i]
			for face in faces:
				t = tuple([j + ii for ii in face])
				outFile.write("3 %d %d %d\n" % t)
			j += len(b.shape.vertices)
		# checks what argument
		self._checkWhatArgumentIsDict(what, "exportPotentialBlocks")
		# write additional data from 'what' param
		if useRef:
			outFile.write("\nPOINT_DATA %d\n" % (len(dspls)))
			outFile.write("\nVECTORS displacement double\n")
			for v in dspls:
				outFile.write("%g %g %g\n" % (v[0], v[1], v[2]))
		if what:
			outFile.write("\nCELL_DATA %d" % (nFaces))
		# see exportSpheres for explanation of this code block
		for name, command in what.items():
			test = eval(command)
			if isinstance(test, Matrix3):
				outFile.write("\nTENSORS %s double\n" % (name))
				for i, b in enumerate(bodies):
					t = eval(command)
					for f in bodyFaces[i]:
						outFile.write(
						        "%g %g %g\n%g %g %g\n%g %g %g\n\n" %
						        (t[0, 0], t[0, 1], t[0, 2], t[1, 0], t[1, 1], t[1, 2], t[2, 0], t[2, 1], t[2, 2])
						)
			elif isinstance(test, Vector3):
				outFile.write("\nVECTORS %s double\n" % (name))
				for i, b in enumerate(bodies):
					v = eval(command)
					for f in bodyFaces[i]:
						outFile.write("%g %g %g\n" % (v[0], v[1], v[2]))
			elif isinstance(test, (int, float)):
				outFile.write("\nSCALARS %s double 1\nLOOKUP_TABLE default\n" % (name))
				for i, b in enumerate(bodies):
					e = eval(command)
					for f in bodyFaces[i]:
						outFile.write("%g\n" % e)
			else:
				self._warn("exportPotentialBlocks: wrong 'what' parameter, vtk output might be corrupted")
		outFile.close()
		self.PotentialBlocksSnapCount += 1


#gmshGeoExport===============================================================
def gmshGeo(filename, comment='', mask=-1, accuracy=-1):
	"""Save spheres in geo-file for the following using in GMSH (http://www.geuz.org/gmsh/doc/texinfo/) program. The spheres can be there meshed.

	:param string filename: the name of the file, where sphere coordinates will be exported.
	:param int mask: export only spheres with the corresponding mask export only spheres with the corresponding mask
	:param float accuracy: the accuracy parameter, which will be set for the poinst in geo-file. By default: 1./10. of the minimal sphere diameter.
	:return: number of spheres which were exported.
	:rtype: int
	"""
	O = Omega()

	try:
		out = open(filename, 'w')
	except:
		raise RuntimeError("Problem to write into the file")

	count = 0
	#out.write('#format \n')
	# Find the minimal diameter
	if (accuracy < 0.0):
		dMin = -1.0
		for b in O.bodies:
			try:
				if (isinstance(b.shape, Sphere) and ((mask < 0) or ((mask & b.mask) > 0))):
					if (((dMin > 0.0) and (dMin > b.shape.radius * 2.0)) or (dMin < 0.0)):
						dMin = b.shape.radius * 2.0
			except AttributeError:
				pass
		accuracy = dMin / 10.0
	# Export bodies
	PTS = 0
	CRS = 0
	out.write('Acc = %g;\n' % (accuracy))
	for b in O.bodies:
		try:
			if (isinstance(b.shape, Sphere) and ((mask < 0) or ((mask & b.mask) > 0))):
				r = b.shape.radius
				x = b.state.pos[0]
				y = b.state.pos[1]
				z = b.state.pos[2]
				out.write('Rad = %g;\n' % (r))
				out.write(
				        'Point(%d) = {%g, %g, %g, Acc};\n\
Point(%d) = {%g, %g, %g, Acc};\n\
Point(%d) = {%g, %g, %g, Acc};\n\
Point(%d) = {%g, %g, %g, Acc};\n\
Point(%d) = {%g, %g, %g, Acc};\n\
Point(%d) = {%g, %g, %g, Acc};\n\
Point(%d) = {%g, %g, %g, Acc};\n\n' % (
				                PTS + 1, x, y, z, PTS + 2, r + x, y, z, PTS + 3, -r + x, y, z, PTS + 4, x, y, r + z, PTS + 5, x, y, -r + z,
				                PTS + 6, x, r + y, z, PTS + 7, x, -r + y, z
				        )
				)
				out.write(
				        '\n\
Circle(%d) = {%d, %d, %d};\n\
Circle(%d) = {%d, %d, %d};\n\
Circle(%d) = {%d, %d, %d};\n\
Circle(%d) = {%d, %d, %d};\n\
Circle(%d) = {%d, %d, %d};\n\
Circle(%d) = {%d, %d, %d};\n\
Circle(%d) = {%d, %d, %d};\n\
Circle(%d) = {%d, %d, %d};\n\
Circle(%d) = {%d, %d, %d};\n\
Circle(%d) = {%d, %d, %d};\n\
Circle(%d) = {%d, %d, %d};\n\
Circle(%d) = {%d, %d, %d};\n' % (
				                CRS + 1,
				                PTS + 4,
				                PTS + 1,
				                PTS + 6,
				                CRS + 2,
				                PTS + 6,
				                PTS + 1,
				                PTS + 5,
				                CRS + 3,
				                PTS + 6,
				                PTS + 1,
				                PTS + 3,
				                CRS + 4,
				                PTS + 3,
				                PTS + 1,
				                PTS + 7,
				                CRS + 5,
				                PTS + 7,
				                PTS + 1,
				                PTS + 5,
				                CRS + 6,
				                PTS + 7,
				                PTS + 1,
				                PTS + 2,
				                CRS + 7,
				                PTS + 2,
				                PTS + 1,
				                PTS + 6,
				                CRS + 8,
				                PTS + 7,
				                PTS + 1,
				                PTS + 4,
				                CRS + 9,
				                PTS + 2,
				                PTS + 1,
				                PTS + 5,
				                CRS + 10,
				                PTS + 5,
				                PTS + 1,
				                PTS + 3,
				                CRS + 11,
				                PTS + 3,
				                PTS + 1,
				                PTS + 4,
				                CRS + 12,
				                PTS + 4,
				                PTS + 1,
				                PTS + 2,
				        )
				)

				out.write(
				        '\n\
Line Loop(%d) = {%d, %d, %d}; Ruled Surface(%d) = {%d};\n\
Line Loop(%d) = {%d, %d, %d}; Ruled Surface(%d) = {%d};\n\
Line Loop(%d) = {%d, %d, %d}; Ruled Surface(%d) = {%d};\n\
Line Loop(%d) = {%d, %d, %d}; Ruled Surface(%d) = {%d};\n\
Line Loop(%d) = {%d, %d, %d}; Ruled Surface(%d) = {%d};\n\
Line Loop(%d) = {%d, %d, %d}; Ruled Surface(%d) = {%d};\n\
Line Loop(%d) = {%d, %d, %d}; Ruled Surface(%d) = {%d};\n\
Line Loop(%d) = {%d, %d, %d}; Ruled Surface(%d) = {%d};\n\n\
' % (
				                (CRS + 13),
				                +(CRS + 1),
				                -(CRS + 7),
				                -(CRS + 12),
				                (CRS + 14),
				                (CRS + 13),
				                (CRS + 15),
				                +(CRS + 7),
				                +(CRS + 2),
				                -(CRS + 9),
				                (CRS + 16),
				                (CRS + 15),
				                (CRS + 17),
				                +(CRS + 2),
				                +(CRS + 10),
				                -(CRS + 3),
				                (CRS + 18),
				                (CRS + 17),
				                (CRS + 19),
				                +(CRS + 3),
				                +(CRS + 11),
				                +(CRS + 1),
				                (CRS + 20),
				                (CRS + 19),
				                (CRS + 21),
				                +(CRS + 8),
				                +(CRS + 12),
				                -(CRS + 6),
				                (CRS + 22),
				                (CRS + 21),
				                (CRS + 23),
				                +(CRS + 4),
				                +(CRS + 8),
				                -(CRS + 11),
				                (CRS + 24),
				                (CRS + 23),
				                (CRS + 25),
				                +(CRS + 5),
				                +(CRS + 10),
				                (CRS + 4),
				                (CRS + 26),
				                (CRS + 25),
				                (CRS + 27),
				                +(CRS + 6),
				                +(CRS + 9),
				                -(CRS + 5),
				                (CRS + 28),
				                (CRS + 27),
				        )
				)
				PTS += 7
				CRS += 28

				count += 1
		except AttributeError:
			pass
	out.close()
	return count


# external vtk manipulation ===============================================================
def text2vtk(inFileName, outFileName, comment="comment"):
	"""Converts text file (created by :yref:`yade.export.textExt` function) into vtk file.
	See :ysrc:`examples/test/paraview-spheres-solid-section/export_text.py` example

	:param str inFileName: name of input text file
	:param str outFileName: name of output vtk file
	:param str comment: optional comment in vtk file
	"""
	with open(inFileName) as fin:
		lines = [line for line in fin]
	for i, line in enumerate(lines):
		if not line.startswith("#"):
			break
		lineColumnNames = line
	columns = lineColumnNames[1:].strip().split()[4:]
	lines = lines[i:]
	data = [line.split() for line in lines]
	n = len(data)
	with open(outFileName, "w") as fout:
		fout.write('# vtk DataFile Version 3.0.\n{}\nASCII\n\nDATASET POLYDATA\nPOINTS {} double\n'.format(comment, n))
		fout.writelines('{} {} {}\n'.format(*d) for d in data)
		fout.write("\nPOINT_DATA {}\nSCALARS radius double 1\nLOOKUP_TABLE default\n".format(n))
		fout.writelines('{}\n'.format(d[3]) for d in data)
		for i, c in enumerate(columns):
			fout.write("\nSCALARS {} double 1\nLOOKUP_TABLE default\n".format(c))
			fout.writelines('{}\n'.format(d[4 + i]) for d in data)


def text2vtkSection(inFileName, outFileName, point, normal=(1, 0, 0)):
	"""Converts section through spheres from text file (created by :yref:`yade.export.textExt` function) into vtk file.
	See :ysrc:`examples/test/paraview-spheres-solid-section/export_text.py` example

	:param str inFileName: name of input text file
	:param str outFileName: name of output vtk file
	:param Vector3|(float,float,float) point: coordinates of a point lying on the section plane
	:param Vector3|(float,float,float) normal: normal vector of the section plane
	"""
	from math import sqrt
	norm = sqrt(pow(normal[0], 2) + pow(normal[1], 2) + pow(normal[2], 2))
	normal = (normal[0] / norm, normal[1] / norm, normal[2] / norm)

	#
	def computeD(point, normal):
		# from point and normal computes parameter d in plane equation ax+by+cz+d=0
		return -normal[0] * point[0] - normal[1] * point[1] - normal[2] * point[2]

	def computeDistanceFromPlane(dat, point, normal, d=None):
		# computes distance of sphere dat from plane (point,normal)
		x, y, z = computeProjectionOnPlane(dat, point, normal, d)
		cx, cy, cz = dat[0], dat[1], dat[2]
		return sqrt(pow(x - cx, 2) + pow(y - cy, 2) + pow(z - cz, 2))

	def computeProjectionOnPlane(self, point, normal, d=None):
		# computes projection of sphere dat on plane (point,normal)
		if d is None:
			d = computeD(point, normal)
		nx, ny, nz = normal[0], normal[1], normal[2]
		cx, cy, cz = dat[0], dat[1], dat[2]
		t = (-d - nx * cx - ny * cy - nz * cz) / (nx * nx + ny * ny + nz * nz)
		x, y, z = cx + t * nx, cy + t * ny, cz + t * nz
		return x, y, z

	#
	fin = open(inFileName)
	lastLine = None
	line = '#'
	while line.startswith('#'):
		lastLine = line
		line = fin.readline()
	columns = lastLine.split()[4:]
	data = [[float(w) for w in line.split()] for line in fin]
	fin.close()
	#
	d = computeD(point, normal)
	circs = []
	for dat in data:
		r = dat[3]
		dst = computeDistanceFromPlane(dat, point, normal, d)
		if dst > r:
			continue
		x, y, z = computeProjectionOnPlane(dat, point, normal, d)
		rNew = sqrt(r * r - dst * dst)
		dNew = [x, y, z, rNew, r]
		dNew.extend(dat[4:])
		circs.append(dNew)
	n = len(circs)
	fout = open(outFileName, 'w')
	fout.write('# vtk DataFile Version 3.0.\ncomment\nASCII\n\nDATASET POLYDATA\nPOINTS %d double\n' % (n))
	fout.writelines('%g %g %g\n' % (c[0], c[1], c[2]) for c in circs)
	fout.write("\nPOINT_DATA %d\nSCALARS radius double 1\nLOOKUP_TABLE default\n" % (n))
	fout.writelines('%g\n' % (c[3]) for c in circs)
	fout.write("\nSCALARS radiusOrig double 1\nLOOKUP_TABLE default\n")
	fout.writelines('%g\n' % (c[4]) for c in circs)
	fout.write("\nVECTORS normal double\n")
	fout.writelines("%g %g %g\n" % normal for i in circs)
	for i, c in enumerate(columns):
		fout.write("\nSCALARS %s double 1\nLOOKUP_TABLE default\n" % (c))
		fout.writelines('%s\n' % (c[4 + i]) for c in circs)
	fout.close()
