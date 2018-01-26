from astropy.table import Table
import glob
from astropy.io.votable import parse_single_table
import matplotlib.pyplot as plt
import pyfits
from sklearn.neighbors import KernelDensity
from sklearn.linear_model import Ridge
from sklearn.linear_model import Lasso
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import GradientBoostingRegressor
from scipy.stats import gaussian_kde
import sqlite3 as lite
import math
import numpy as np
import seaborn as sns
from pandas.tools.plotting import scatter_matrix
import pandas as pd
from sklearn.model_selection import KFold
sns.set(font_scale=4)


def createandfilltables():
	con = lite.connect('/home/hartsuiker/Documents/dbdm/DDM2017/FinalProject/DDM17final.db')

	with con:
		cur = con.cursor()

		for fitsName in glob.glob('/home/hartsuiker/Documents/dbdm/DDM2017/FinalProject/Q1/Tables/*.csv'):
			command = """CREATE TABLE IF NOT EXISTS {0} (ID INT,
			FieldID INT,Filename varchar(30),Filter varchar(5), MJD DOUBLE, 
			Airmass FLOAT, Exptime FLOAT,UNIQUE(ID), PRIMARY KEY(ID))""".format('mastertable')
			cur.execute(command)
			con.commit()
			fitsdata = Table.read(fitsName)
			for row in fitsdata:
				command2="INSERT OR IGNORE INTO mastertable VALUES({0},{1},'{2}','{3}',{4},{5},{6}) ".format(row[0],row[1],row[2],row[3],row[4],row[5],row[6])
				cur.execute(command2)
			con.commit()

		AAA = np.vstack((fitsdata[2][:],fitsdata[0][:],fitsdata[3][:]))
		BBB = list(set(AAA[3][:]))

		cur.execute(command)
		for i in range(len(BBB)):
			command = """CREATE TABLE IF NOT EXISTS {0} (ImageID INT,RunningID FLOAT,
			X FLOAT , Y FLOAT, Flux1 FLOAT,dFlux1 FLOAT, Flux2 FLOAT,dFlux2 FLOAT,
			Flux3 FLOAT, dFlux3 FLOAT, Ra FLOAT, Dec FLOAT, Class INT,Mag1 FLOAT,dMag1 FLOAT
			,Mag2 FLOAT,dMag2 FLOAT,Mag3 FLOAT,dMag3 FLOAT,StarID INT)""".format('imagetable_'+str(BBB[i]))
			cur.execute(command)
			con.commit()
			for fitsName in glob.glob('/home/hartsuiker/Documents/dbdm/DDM2017/FinalProject/Q1/Tables/Field-*'+str(BBB[i])+'*.fits'):
				for k in range(len(AAA[0])):
					if fitsName == '/home/hartsuiker/Documents/dbdm/DDM2017/FinalProject/Q1/Tables/'+AAA[0][k]:
						imageid = AAA[1][k]
				fitsdata = Table.read(fitsName)
				for row in fitsdata:
					command3="INSERT INTO imagetable_"+str(BBB[i])+" VALUES("+str(imageid)+","
					for j in range(len(fitsdata[0])):
						if math.isnan(row[j]):
							command3 += "NULL"
						else:
							command3 += str(row[j])
						if j != (len(fitsdata[0])-1):
							command3 += ","
					command3 += ")"
					cur.execute(command3)
			con.commit()

def doqueries(givenfield,command):
	con = lite.connect('/home/hartsuiker/Documents/dbdm/DDM2017/FinalProject/DDM17final.db')

	with con:
		cur = con.cursor()

		commandR1 = """SELECT ImageID,COUNT(DISTINCT StarID)
						FROM imagetable_H
						WHERE Flux1/dFlux1 > 5
						and ImageID in(
							SELECT ID
							FROM mastertable
							WHERE MJD between 56800 and 57300)
						and class = -1
						GROUP BY ImageID

						UNION

						SELECT ImageID,COUNT(DISTINCT StarID)
						FROM imagetable_Ks
						WHERE Flux1/dFlux1 > 5
						and ImageID in(
							SELECT ID
							FROM mastertable
							WHERE MJD between 56800 and 57300)
						and class = -1
						GROUP BY ImageID

						UNION

						SELECT ImageID,COUNT(DISTINCT StarID)
						FROM imagetable_Z
						WHERE Flux1/dFlux1 > 5
						and ImageID in(
							SELECT ID
							FROM mastertable
							WHERE MJD between 56800 and 57300)
						and class = -1
						GROUP BY ImageID

						UNION

						SELECT ImageID,COUNT(DISTINCT StarID)
						FROM imagetable_J
						WHERE Flux1/dFlux1 > 5
						and ImageID in(
							SELECT ID
							FROM mastertable
							WHERE MJD between 56800 and 57300)
						and class = -1
						GROUP BY ImageID

						UNION

						SELECT ImageID,COUNT(DISTINCT StarID)
						FROM imagetable_Y
						WHERE Flux1/dFlux1 > 5
						and ImageID in(
							SELECT ID
							FROM mastertable
							WHERE MJD between 56800 and 57300)
						and class = -1
						GROUP BY ImageID
						ORDER BY imageID asc
						"""

		commandR2 = '''SELECT h.StarID,j.mag1-h.mag1
						FROM imagetable_H as h
						join imagetable_J as j on h.StarID = j.StarID
						WHERE J.mag1-h.mag1 > 1.5
						ORDER BY h.StarID asc
					'''

		commandR3 = '''SELECT ks.StarID,ks.imageID,ABS(ks.Flux1-(
														SELECT AVG(ks2.Flux1)
														FROM imagetable_Ks as ks2
														WHERE ks.imageID = ks2.imageID)
														)/ks.dFlux1
						FROM imagetable_Ks as ks
						WHERE ABS(ks.Flux1 -(
							SElECT AVG(ks2.Flux1)
							FROM imagetable_Ks as ks2
							WHERE ks.imageID = ks2.imageID)) > 20 *ks.dFlux1
						ORDER BY ks.StarID asc,ks.imageID asc
					'''

		commandR4 = '''SELECT ID
						FROM mastertable
						WHERE FieldID = %s
						ORDER BY ID asc
					'''%(givenfield)

		commandR5 = '''SELECT y.StarID,y.Mag1,z.Mag1,j.Mag1,h.Mag1,ks.Mag1
						FROM imagetable_Y as y
						join imagetable_Z as z on z.StarID = y.StarID
						join imagetable_J as j on j.StarID = y.StarID
						join imagetable_H as h on h.StarID = y.StarID
						join imagetable_Ks as ks on ks.StarID = y.StarID
						join mastertable as m on m.ID = y.ImageID
						WHERE y.Flux1/y.dFlux1 > 30
						and z.Flux1/z.dFlux1 > 30
						and j.Flux1/j.dFlux1 > 30
						and h.Flux1/h.dFlux1 > 30
						and ks.Flux1/ks.dFlux1 > 30
						and y.class = -1
						and z.class = -1
						and j.class = -1
						and h.class = -1
						and ks.class = -1
						and m.FieldID = %s
						and ks.ImageID=(
							SELECT m2.ID
							FROM mastertable as m2
							WHERE m2.Filename = 'Field-%s-Ks-E001.fits')
						ORDER BY y.StarID asc
					'''%(givenfield,givenfield)

		commandR6 = '''SELECT y.Mag1-j.Mag1,j.Mag1-h.Mag1
						FROM imagetable_Y as y
						join imagetable_J as j on j.StarID = y.StarID
						join imagetable_H as h on h.StarID = y.StarID
						WHERE y.Mag1-j.Mag1 not NULL
						and j.Mag1-h.Mag1 not NULL
						and y.class = -1
						and j.class = -1
						and h.class = -1
						limit 100
					'''
		if command ==1:
			command = commandR1
		elif command == 2:
			command = commandR2
		elif command == 3:
			command = commandR3
		elif command == 4:
			command = commandR4
		elif command == 5:
			command = commandR5
		elif command == 6:
			command = commandR6
		rows=cur.execute(command)

		# for row in rows:
		# 	print row

		if command == commandR2:
			Q=0
			for row in rows:
				if Q==0:
					a=np.array(row)
					Q=1
				else:
					a=np.vstack((a,row))
			print a
			plt.hist(a[:,1],bins=200)
			plt.ylabel('amount of objects',fontsize=50)
			plt.xlabel('J-H color',fontsize=50)
			plt.xticks(fontsize=40)
			plt.yticks(fontsize=40)
			plt.xlim(1.49,1.75)
			plt.title('J-H color of all objects with J-H > 1.5',fontsize=60)
			plt.show()
			plt.close()

		if command == commandR3:
			Q=0
			for row in rows:
				if Q==0:
					a=np.array(row)
					Q=1
				else:
					a=np.vstack((a,row))
			print a
			plt.hist(a[:,2],bins=280)
			plt.ylabel('amount of objects',fontsize=50)
			plt.xlabel('deviation from the mean flux [flux uncertainties]',fontsize=50)
			plt.xticks(fontsize=40)
			plt.yticks(fontsize=40)
			plt.xlim(0,145)
			plt.title('deviation from the mean flux for all deviations > 20 times the flux uncertainty',fontsize=32)
			plt.show()
			plt.close()

		if command == commandR5:
			Q=0
			for row in rows:
				if Q==0:
					a=np.array(row)
					Q=1
				else:
					a=np.vstack((a,row))

			sns.kdeplot(a[:,1],label='Y',shade=True,linewidth=3.5)
			sns.kdeplot(a[:,2],label='Z',shade=True,linewidth=3.5)
			sns.kdeplot(a[:,3],label='J',shade=True,linewidth=3.5)
			sns.kdeplot(a[:,4],label='H',shade=True,linewidth=3.5)
			sns.kdeplot(a[:,5],label='Ks',shade=True,linewidth=3.5)
			plt.title('Kernel density plot in all filters of all objects in field '+str(givenfield),fontsize=50)
			leg = plt.legend(fontsize=60,loc='upper left')
			for line in leg.get_lines():
				line.set_linewidth(6.0)
			plt.xlabel('Magnitude in given filter',fontsize=50)
			plt.ylabel('Normalized counts',fontsize=50)
			plt.xticks(fontsize=30)
			plt.yticks(fontsize=30)
			plt.show()
			plt.close()


		if command == commandR6:
			Q=0
			for row in rows:
				if Q==0:
					a=np.array(row)
					Q=1
				else:
					a=np.vstack((a,row))


			kf = KFold(n_splits=10)
			kf.get_n_splits(a)
			print 'shape',np.shape(a)
			Max=-1e99
			# for i in range(1000):
			# 	print 0.001+i/1000.
			# 	array=[]
			# 	for train_index,test_index in kf.split(a):
			# 		a_train,a_test = a[train_index],a[test_index]
			# 		kde = KernelDensity(kernel='gaussian', bandwidth=0.001+i/1000.).fit(a_train)
			# 		log_dens = kde.score_samples(a_train)
			# 		loglikelihood = kde.score(a_test)
			# 		array = np.append(array,loglikelihood)
			# 	Loglikelihood = np.nanmean(array)
			# 	if Loglikelihood > Max:
			# 		Max=Loglikelihood
			# 		Bandwidth = 0.001+i/1000.
			# 		print 'new best value for the bandwidth: ',Bandwidth
			Bandwidth=0.061 #calculated with the above for loop for the fist 2000 entries of the query
			kde = KernelDensity(kernel='gaussian', bandwidth=Bandwidth).fit(a)
			samples = kde.sample(100000)
			# plt.scatter(samples[:,0],samples[:,1])
			# plt.xlabel('Y-J',fontsize=50)
			# plt.ylabel('J-H',fontsize=50)
			# plt.xticks(fontsize=40)
			# plt.yticks(fontsize=40)
			# plt.title('sample of J-H color vs the Y-J color for 100,000 stars',fontsize=50)
			# plt.show()
			# plt.close()
			data = samples
			df = pd.DataFrame(data, columns=["Y-J", "J-H"])
			sns.jointplot(x="Y-J", y="J-H", data=df, stat_func = None, kind="kde")
			# plt.xlabel('Y-J',fontsize=40)
			# plt.ylabel('J-H',fontsize=40)
			plt.xticks(fontsize=20)
			plt.yticks(fontsize=20)
			# plt.title('sample of J-H color vs the Y-J color for 100,000 stars as 2D distribution',fontsize=40)
			plt.show()
			plt.close()
		con.commit()

# createandfilltables()

givenfield = 3
doqueries(givenfield,6)
kjhkjh

tableA = parse_single_table('/home/hartsuiker/Documents/dbdm/DDM2017/FinalProject/Q2/Tables/PhotoZFileA.vot')
arrayA = np.array([tableA.array['mag_r'],tableA.array['u-g'],tableA.array['g-r'],tableA.array['r-i'],tableA.array['i-z'],tableA.array['z_spec']]).T

g = sns.pairplot(pd.DataFrame(arrayA[:500,0:5],columns=('mag_r','u-g','g-r','r-i','i-z')))
plt.show()
plt.close()

def linearregressiontrain():
	regr = LinearRegression()

	result = regr.fit(arrayA[:,0:-1],arrayA[:,-1])
	Zpredtrain = regr.predict(arrayA[:,0:-1])
	E = np.median(abs((arrayA[:,-1]-Zpredtrain)/(1+arrayA[:,-1])))
	return E

def ridgeregressiontrain(Alpha):
	regr = Ridge(alpha=Alpha)

	result = regr.fit(arrayA[:,0:-1],arrayA[:,-1])
	Zpredtrain = regr.predict(arrayA[:,0:-1])
	E = np.median(abs((arrayA[:,-1]-Zpredtrain)/(1+arrayA[:,-1])))
	return E

def lassoregressiontrain(Alpha):
	regr = Lasso(alpha=Alpha)

	result = regr.fit(arrayA[:,0:-1],arrayA[:,-1])
	Zpredtrain = regr.predict(arrayA[:,0:-1])
	E = np.median(abs((arrayA[:,-1]-Zpredtrain)/(1+arrayA[:,-1])))
	return E

def linearregressiontrain():
	Ermin=Elamin=1e99
	alphaEr=alphaEla=np.pi
	for i in range(1,50000):
		alpha = i/100000.
		print alpha
		Er = ridgeregressiontrain(alpha)
		Ela = lassoregressiontrain(alpha)
		if Er < Ermin:
			Ermin=Er
			alphaEr = alpha
		if Ela < Elamin:
			Elamin=Ela
			alphaEla = alpha
	Elimin = linearregressiontrain()
	print 'linear regression training error: ',Elimin
	print 'ridge regression training error: ',Ermin,' at alpha = ',alphaEr
	print 'lasso regression training error: ',Elamin,' at alpha = ',alphaEla

tableB = parse_single_table('/home/hartsuiker/Documents/dbdm/DDM2017/FinalProject/Q2/Tables/PhotoZFileB.vot')
arrayB = np.array([tableB.array['mag_r'],tableB.array['u-g'],tableB.array['g-r'],tableB.array['r-i'],tableB.array['i-z'],tableB.array['z_spec']]).T

def linearregressiontest():	
	regr2 = LinearRegression()
	result = regr2.fit(arrayA[:,0:-1],arrayA[:,-1])
	Zpredtest = regr2.predict(arrayB[:,0:-1])

	E = np.median(abs((arrayB[:,-1]-Zpredtest)/(1+arrayB[:,-1])))
	print 'linear test ',E

def lassoregressiontest():	
	regr2 = Lasso(alpha =8e-7)
	result = regr2.fit(arrayA[:,0:-1],arrayA[:,-1])
	Zpredtest = regr2.predict(arrayB[:,0:-1])

	E = np.median(abs((arrayB[:,-1]-Zpredtest)/(1+arrayB[:,-1])))
	print 'lasso test ',E

def ridgeregressiontest():	
	regr2 = Ridge(alpha =0.1905)
	result = regr2.fit(arrayA[:,0:-1],arrayA[:,-1])
	Zpredtest = regr2.predict(arrayB[:,0:-1])

	E = np.median(abs((arrayB[:,-1]-Zpredtest)/(1+arrayB[:,-1])))
	print 'ridge test: ',E

linearregressiontest()
lassoregressiontest()
ridgeregressiontest()


def bestbandwidth(a):
	kf = KFold(n_splits=10)
	kf.get_n_splits(a)
	Max=-1e99
	for train_index,test_index in kf.split(a):
		a_train,a_test = a[train_index],a[test_index]
		kde = KernelDensity(kernel='gaussian', bandwidth=0.001+i/1000.).fit(a_train)
		log_dens = kde.score_samples(a_train)
		loglikelihood = kde.score(a_test)
		array = np.append(array,loglikelihood)
	Loglikelihood = np.nanmean(array)
	if Loglikelihood > Max:
		Max=Loglikelihood
		Bandwidth = 0.001+i/1000.
		print 'new best value for the bandwidth: ',Bandwidth
bestbandwidth(a)



def kneighbours():
	Emin = 1e90
	imin = 0
	numberofsplits = 10
	kf = KFold(n_splits=numberofsplits)
	kf.get_n_splits(arrayA)
	for i in range(1,30):
		E=0
		for train_index,test_index in kf.split(arrayA):
			a_train,a_test = arrayA[train_index],arrayA[test_index]
			neigh = KNeighborsRegressor(n_neighbors=i)
			result = neigh.fit(a_train[:,0:-1],a_train[:,-1])
			Zpredtest = neigh.predict(a_test[:,0:-1])
			E += np.median(abs((a_test[:,-1]-Zpredtest)/(1+a_test[:,-1])))/numberofsplits
		if E < Emin:
			Emin = E
			imin = i
	imin=19
	print 'minimal error on the training set is ',Emin,'for ',imin,'nearest neighbours'
	neigh = KNeighborsRegressor(n_neighbors=19)
	result = neigh.fit(arrayA[:,0:-1],arrayA[:,-1]) 
	Zpred = neigh.predict(arrayB[:,0:-1])	
	E = np.median(abs((arrayB[:,-1]-Zpred)/(1+arrayB[:,-1])))
	print 'the generalized error on the test set is: ',E

kneighbours()

def randomforest():
	Emin = 1e90
	imin = 0
	numberofsplits = 10
	kf = KFold(n_splits=numberofsplits)
	kf.get_n_splits(arrayA)
	# for i in range(1,10):
	# 	E=0
	# 	for train_index,test_index in kf.split(arrayA):
	# 		a_train,a_test = arrayA[train_index],arrayA[test_index]
	# 		neigh = RandomForestRegressor(n_estimators=100+10*i,n_jobs=-3)
	# 		result = neigh.fit(a_train[:,0:-1],a_train[:,-1])
	# 		Zpredtest = neigh.predict(a_test[:,0:-1])
	# 		E += np.median(abs((a_test[:,-1]-Zpredtest)/(1+a_test[:,-1])))/numberofsplits
	# 	if E < Emin:
	# 		Emin = E
	# 		imin = i
	imin=130
	print 'minimal error on the training set is ',Emin,'for ',imin,' trees'
	neigh = RandomForestRegressor(n_estimators=imin)
	result = neigh.fit(arrayA[:,0:-1],arrayA[:,-1]) 
	Zpred = neigh.predict(arrayB[:,0:-1])	
	E = np.median(abs((arrayB[:,-1]-Zpred)/(1+arrayB[:,-1])))
	print 'the generalized error on the test set is: ',E

randomforest()



def neuralnetwork():
	Emin = 1e90
	imin = 0
	numberofsplits = 10
	kf = KFold(n_splits=numberofsplits)
	kf.get_n_splits(arrayA)
	for i in range(1,30):
		E=0
		for train_index,test_index in kf.split(arrayA):
			a_train,a_test = arrayA[train_index],arrayA[test_index]
			neigh = MLPRegressor(hidden_layer_sizes=(i*10, ))
			result = neigh.fit(a_train[:,0:-1],a_train[:,-1])
			Zpredtest = neigh.predict(a_test[:,0:-1])
			E += np.median(abs((a_test[:,-1]-Zpredtest)/(1+a_test[:,-1])))/numberofsplits
		if E < Emin:
			Emin = E
			imin = i*10
	print 'minimal error on the training set is ',Emin,'for a hidden layer size of',imin
	neuraln = MLPRegressor(hidden_layer_sizes=(imin, ))
	result = neuraln.fit(arrayA[:,0:-1],arrayA[:,-1]) 
	Zpred = neuraln.predict(arrayB[:,0:-1])
	E = np.median(abs((arrayB[:,-1]-Zpred)/(1+arrayB[:,-1])))
	print 'neural network: ',E

neuralnetwork()

def gradientboosting():
	Emin = 1e90
	imin = 0
	numberofsplits = 10
	kf = KFold(n_splits=numberofsplits)
	kf.get_n_splits(arrayA)
	for i in range(1,70):
		E=0
		for train_index,test_index in kf.split(arrayA):
			a_train,a_test = arrayA[train_index],arrayA[test_index]
			neigh = GradientBoostingRegressor(n_estimators=5*i)
			result = neigh.fit(a_train[:,0:-1],a_train[:,-1])
			Zpredtest = neigh.predict(a_test[:,0:-1])
			E += np.median(abs((a_test[:,-1]-Zpredtest)/(1+a_test[:,-1])))/numberofsplits
		if E < Emin:
			Emin = E
			imin = 5*i
			print Emin,imin
	print 'minimal error on the training set is ',Emin,'for ',imin,' boosting stages'
	grboost = GradientBoostingRegressor(n_estimators = imin)
	result = grboost.fit(arrayA[:,0:-1],arrayA[:,-1]) 
	Zpred = grboost.predict(arrayB[:,0:-1])
	E = np.median(abs((arrayB[:,-1]-Zpred)/(1+arrayB[:,-1])))
	print 'gradient boosting regressor: ',E

gradientboosting()
