# -*- encoding=utf-8 -*-
# coding: utf-8
# 2009 © Václav Šmilauer <eudoxos@arcig.cz>
# 2011 ©Bruno Chareyre <bruno.chareyre@grenoble-inp.fr>
"Test and demonstrate use of PeriTriaxController."
from yade import pack, qt

O.periodic = True
O.cell.setBox(.1, .1, .1)
#O.cell.Hsize=Matrix3(0.1,0,0, 0,0.1,0, 0,0,0.1)
sp = pack.SpherePack()
radius = 5e-3
num = sp.makeCloud((0, 0, 0), (.1, .1, .1), radius, .2, 500, periodic=True)  # min,max,radius,rRelFuzz,spheresInCell,periodic
O.bodies.append([sphere(s[0], s[1]) for s in sp])

O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb()], verletDist=.05 * radius),
        InteractionLoop([Ig2_Sphere_Sphere_ScGeom()], [Ip2_FrictMat_FrictMat_FrictPhys()], [Law2_ScGeom_FrictPhys_CundallStrack()]),
        PeriTriaxController(
                dynCell=True,
                mass=0.2,
                maxUnbalanced=0.01,
                relStressTol=0.02,
                goal=[-1e4, -1e4, 0],
                stressMask=3,
                globUpdate=5,
                maxStrainRate=[1., 1., 1.],
                doneHook='triaxDone()',
                label='triax'
        ),
        NewtonIntegrator(damping=.2),
]

phase = 0


def triaxDone():
	global phase
	if phase == 0:
		print('Here we are: stress', triax.stress, 'strain', triax.strain, 'stiffness', triax.stiff)
		print('\nNow shearing. Press ▶ (the start button) to start shearing.\n')
		O.cell.velGrad = Matrix3(0, 6, 0, 0, 0, 0, 0, 0, 0)
		triax.stressMask = 7
		triax.goal = [-1e4, -1e4, -1e4]
		phase += 1
		O.pause()
	#elif phase==1:
	#print 'Here we are: stress',triax.stress,'strain',triax.strain,'stiffness',triax.stiff
	#phase+=1
	##print 'Done, pausing now.'
	##O.pause()


O.dt = PWaveTimeStep()

#O.run(7000);
#qt.View()
##r=qt.Renderer()
##r.bgColor=1,1,1
#O.wait()
#O.saveTmp()
#O.cell.velGrad=Matrix3(0,0,0, -6,0,0, 0,0,0)
#O.run(5000);
#O.wait()

#O.cell.velGrad=Matrix3(0,6,0, 0,0,0, 0,0,0)
#O.run(5000);
