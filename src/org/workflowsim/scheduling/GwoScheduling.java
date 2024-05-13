package org.workflowsim.scheduling;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.Map;

public class GwoScheduling {
	
//	public static int taskNum = 3;
//	public static int particleNum = 10;
//	public static int vmNum = 3;
//	public static int iterateNum = 10000;
//	
//	public static double a = 2.0;
//	public static double[] r1 = new double[particleNum];
//	public static double[] r2 = new double[particleNum];
//	public static double[] A = new double[particleNum];
//	public static double[] C = new double[particleNum];
//	
//	public static double[] X_alpha = new double[taskNum];
//	public static double[] X_beta = new double[taskNum];
//	public static double[] X_delta = new double[taskNum];
//	public static List<double[]> velocity=new ArrayList<double[]>();
//	public static double alpha_best_fitness=Double.MAX_VALUE;
//	public static double beta_best_fitness=Double.MAX_VALUE;
//	public static double delta_best_fitness=Double.MAX_VALUE;
//	public static int[] x;
//	public static int[] pbestSchedule;
//	public static int initFlag=0;
//	public static List<int[]> schedules=new ArrayList<int[]>();
//	public static List<int[]> newSchedules=new ArrayList<int[]>();
//	public static List<int[]> alpha_best_schedule=new ArrayList<int[]>();
//	public static List<int[]> beta_best_schedule=new ArrayList<int[]>();
//	public static List<int[]> delta_best_schedule=new ArrayList<int[]>();
//	public static int[] alphabestSchedule;
//	public static int[] betabestSchedule;
//	public static int[] deltabestSchedule;
//	public static double[] pbest_fitness;
//	public static List<int[]> pbest_schedule=new ArrayList<int[]>();
////	public static double[] execcost = new double[vmNum];
////	public static double[] waitcost = new double[taskNum];
////	public static int[] outputfilesize = new int[taskNum];
////	public static int[] tasklength = new int[taskNum];
////	public static int[] mips = new int[vmNum];
////	public static int[][] graph = new int[taskNum][taskNum];
//	
//	public static void init(int jobNum,int maxVmNum) {
//		taskNum = jobNum;
//		vmNum = maxVmNum;
////		initializeParticles();
//		initializeValues();
////		initializeMatrices();
//		pbest_fitness = new double[particleNum];
//		alphabestSchedule = new int[taskNum];
//		betabestSchedule = new int[taskNum];
//		deltabestSchedule = new int[taskNum];
//		Random random = new Random();
//		
//		for(int i=0;i<particleNum;i++)
//		{
//			x=new int[taskNum];
//			pbestSchedule=new int[taskNum];
//
//			for(int j=0;j<taskNum;j++)
//			{
//				x[j] = (random.nextInt(vmNum) % vmNum);
//				pbestSchedule[j]=x[j];
//			}
//			schedules.add(x);
//			pbest_schedule.add(pbestSchedule);
//		}
//		initFlag=1;
////		Random random = new Random();
////       
////		FitnessFunctionGWO fitness = new FitnessFunctionGWO(execcost,waitcost,mips,outputfilesize,tasklength,graph,taskNum,vmNum);
////		updateBestParticles(fitness);
////		int t = 0;
////		while(t < iterateNum) {
////			for(int i = 0; i < particleNum; i++) {
////				particles[i] = updatePosition(i);
//// 			}
////			updateValues();
////			updateBestParticles(fitness);
////			t++;
////		}
//	}
//	
////	public static void initializeParticles() {
////		Random random = new Random();
////		for(int i = 0; i < particleNum; i++) {
////			for(int j = 0; j < taskNum; j++) {
////				particles[i][j] = (random.nextInt(vmNum) % vmNum);				
////			}
////		}
////	}
//	
//	public static void initializeValues() {
//		Random random = new Random();
//		for(int i = 0; i < particleNum; i ++) {
//			double randomDouble = random.nextDouble();
//			r1[i] = randomDouble;
//			randomDouble = random.nextDouble();
//			r2[i] = randomDouble;
//			A[i] = 2*a*r1[i] - a;
//			C[i] = 2*r2[i];
//		}
//	}
//	
////	public static void initializeMatrices() {
////		Random random = new Random();
////		for(int i = 0; i < vmNum; i++) {
////			execcost[i] = random.nextDouble() * 100;
////			mips[i] = (int)Math.round(random.nextDouble() * 100000000);
////		}
////		for(int i = 0; i < taskNum; i++) {
////			waitcost[i] = random.nextDouble() * 100;
////			outputfilesize[i] = (int)Math.round(random.nextDouble() * 10000);
////			tasklength[i] = (int)Math.round(random.nextDouble() * 1000);
////			for(int j = 0; j < taskNum; j++) {
////				graph[i][j] = 0;
////			}
////		}
////	}
////	
////	public static void updateBestParticles(FitnessFunctionGWO fitness) {
////		double[] currentCost = new double[particleNum];
////		for(int i = 0; i < particleNum; i++) {
////			currentCost[i] = fitness.calculatecost(particles[i]);
////		}
////		
////		int[] minKeys = new int[3];
////		
////		minKeys[0] = -1;
////		minKeys[1] = -1;
////		minKeys[2] = -1;
////		
////        for (int i = 0; i < 3; i++) {
////            double minValue = Double.MAX_VALUE;
////            int minKey = 0;
////            for (int j = 0; j < particleNum; j++) {
////                int key = j;
////                double value = currentCost[j];
////                if (value < minValue && minKeys[0] != key && minKeys[1] != key && minKeys[2] != key) {
////                    minValue = value;
////                    minKey = key;
////                }
////            }
////            minKeys[i] = minKey;
////        }
////        
////        X_alpha = particles[minKeys[0]];
////        X_beta = particles[minKeys[1]];
////        X_delta = particles[minKeys[2]];
////	}
//	
//	public static void updateValues()
//	{
//		a -= (2 / iterateNum);
//		Random random = new Random();
//		for(int i = 0; i < particleNum; i ++) {
//			double randomDouble = random.nextDouble();
//			r1[i] = randomDouble;
//			randomDouble = random.nextDouble();
//			r2[i] = randomDouble;
//			A[i] = 2*a*r1[i] - a;
//			C[i] = 2*r2[i];
//		}
//	}
//	
//	public static void updateParticles() {
//		Random random = new Random();
//		for(int i = 0; i < particleNum; i++) {
//			int x[]=schedules.get(i);
//			int alphabest[]=new int[taskNum];
//			int betabest[]=new int[taskNum];
//			int deltabest[]=new int[taskNum];
//			
//			alphabest=alphabestSchedule;
//			betabest=betabestSchedule;
//			deltabest=deltabestSchedule;
//			
//			for(int j = 0; j < x.length; j++) {
//				double D_alpha = Math.abs(C[random.nextInt(particleNum)] * alphabest[j] - x[j]);
//				double D_beta = Math.abs(C[random.nextInt(particleNum)] * betabest[j] - x[j]);
//				double D_delta = Math.abs(C[random.nextInt(particleNum)] * deltabest[j] - x[j]);
//				double X1 = alphabest[j] - A[random.nextInt(particleNum)] * D_alpha;
//				double X2 = betabest[j] - A[random.nextInt(particleNum)] * D_beta;
//				double X3 = deltabest[j] - A[random.nextInt(particleNum)] * D_delta;
//				x[j] = (Math.abs((int)Math.round((X1+X2+X3)/3)))%vmNum;
//			}			
//		}
//		newSchedules.add(x);
//		updateValues();
//	}
//	
//	/**
//	 * Initialize all objects, in order to repeatedly implement the pso adjustment algorithm
//	 */
//	public static void clear() {
//		alpha_best_fitness = Double.MAX_VALUE;
//		beta_best_fitness = Double.MAX_VALUE;
//		delta_best_fitness = Double.MAX_VALUE;
//	    initFlag = 0;
//        schedules.removeAll(schedules);
//        pbest_schedule.removeAll(pbest_schedule);
//        newSchedules.removeAll(newSchedules);
//        pbest_schedule.removeAll(pbest_schedule);
//	}
	
	public static int particleNum = 10000;//粒子数
	public static int iterateNum = 3000;//迭代次数
	public static double c1 = 2.0;//学习因子c1
	public static double c2 = 1.5;//学习因子c2
	public static double w = 3;//惯性权重
	public static double a = 2.0;
	public static double[] r1 = new double[particleNum];
	public static double[] r2 = new double[particleNum];
	public static double[] A = new double[particleNum];
	public static double[] C = new double[particleNum];
	public static int initFlag = 0;
	public static List<int[]> schedules=new ArrayList<int[]>();//更新前
	public static List<int[]> newSchedules=new ArrayList<int[]>();//更新后
	public static List<int[]> pbest_schedule=new ArrayList<int[]>();
	public static int[] alphabest_schedule;
	public static int[] betabest_schedule;
	public static int[] deltabest_schedule;
	public static double[] pbest_fitness;
	public static List<double[]> velocity=new ArrayList<double[]>();
	public static double alphabest_fitness=Double.MAX_VALUE;
	public static double betabest_fitness=Double.MAX_VALUE;
	public static double deltabest_fitness=Double.MAX_VALUE;
	public static int[] x;//更新前
	public static int[] pbestSchedule;
	public static double[] v;
	public static int taskNum;
	public static int vmNum;
	
	public static void init(int jobNum,int maxVmNum) {
		initializeValues();
		pbest_fitness=new double[particleNum];
		  taskNum=jobNum;
		  vmNum=maxVmNum;
		  alphabest_schedule=new int[taskNum];
		  betabest_schedule=new int[taskNum];
		  deltabest_schedule=new int[taskNum];
			for(int i=0;i<particleNum;i++)
			{
				x=new int[taskNum];
				pbestSchedule=new int[taskNum];
				double[] v=new double[taskNum];
				for(int j=0;j<taskNum;j++)
				{
					x[j]=new Random().nextInt(vmNum);
					pbestSchedule[j]=x[j];
					v[j]=new Random().nextDouble();
				}
				schedules.add(x);
				pbest_schedule.add(pbestSchedule);
				velocity.add(v);
			}
			initFlag=1;
		}
	
	public static void initializeValues() {
	Random random = new Random();
	for(int i = 0; i < particleNum; i ++) {
		double randomDouble = random.nextDouble();
		r1[i] = randomDouble;
		randomDouble = random.nextDouble();
		r2[i] = randomDouble;
		A[i] = 2*a*r1[i] - a;
		C[i] = 2*r2[i];
	}
}	
	public static void updateValues()
	{
		a -= (2 / iterateNum);
		Random random = new Random();
		for(int i = 0; i < particleNum; i ++) {
			double randomDouble = random.nextDouble();
			r1[i] = randomDouble;
			randomDouble = random.nextDouble();
			r2[i] = randomDouble;
			A[i] = 2*a*r1[i] - a;
			C[i] = 2*r2[i];
		}
	}
	
//	public static void updateParticles()
//	{
//		for(int i=0;i<particleNum;i++) {
//			int x[]=schedules.get(i);
//			double v[]=velocity.get(i);
//			int pbest[]=new int[taskNum];
//			int temp1[]=new int[taskNum];
//			int temp2[]=new int[taskNum];
//			double sum[]=new double[taskNum];
//			//更新每种调度方案
//			pbest=pbest_schedule.get(i);
//			for(int j=0;j<x.length;j++) {
//				temp1[j] = pbest[j]-x[j];
//				temp2[j] = alphabest_schedule[j]-x[j];
//				double r1 = new Random().nextDouble();
//				double r2 = new Random().nextDouble();
//				sum[j] = c1*r1*temp1[j]+c2*r2*temp2[j];
//				v[j] = w*v[j]+sum[j];
//				x[j] = x[j]+(int)v[j];
//				if(x[j]>vmNum-1)
//					x[j]=vmNum-1;
//				if(x[j]<0)
//					x[j]=0;
//			}
//			newSchedules.add(x);
//		}
//	}
	
	public static void updateParticles() {
	Random random = new Random();
	for(int i = 0; i < particleNum; i++) {
		int x[]=schedules.get(i);
		int alphabest[]=new int[taskNum];
		int betabest[]=new int[taskNum];
		int deltabest[]=new int[taskNum];
		double sum[]=new double[taskNum];
		double v[]=velocity.get(i);
		int pbest[]=new int[taskNum];
		int temp1[]=new int[taskNum];
		int temp2[]=new int[taskNum];
		alphabest=alphabest_schedule;
		betabest=betabest_schedule;
		deltabest=deltabest_schedule;
		pbest=pbest_schedule.get(i);

		for(int j = 0; j < x.length; j++) {
			double D_alpha = Math.abs(C[random.nextInt(particleNum)] * alphabest[j] - x[j]);
			double D_beta = Math.abs(C[random.nextInt(particleNum)] * betabest[j] - x[j]);
			double D_delta = Math.abs(C[random.nextInt(particleNum)] * deltabest[j] - x[j]);
			double X1 = alphabest[j] - A[random.nextInt(particleNum)] * D_alpha;
			double X2 = betabest[j] - A[random.nextInt(particleNum)] * D_beta;
			double X3 = deltabest[j] - A[random.nextInt(particleNum)] * D_delta;
			x[j] = (Math.abs((int)Math.round((X1+X2+X3)/3)))%vmNum;
			if(x[j] >= vmNum) {
				x[j] = vmNum - 1; 
			}
			if(x[j] < 0) {
				x[j] = 0;
			}
		}			
	}
	
	newSchedules.add(x);
	updateValues();
}
	
	/**
	 * 初始化所有对象，为了反复实现pso调度算法
	 */
	public static void clear() {
		alphabest_fitness = Double.MAX_VALUE;
	    initFlag = 0;
        schedules.removeAll(schedules);
        pbest_schedule.removeAll(pbest_schedule);
        velocity.removeAll(velocity);
        newSchedules.removeAll(newSchedules);
        pbest_schedule.removeAll(pbest_schedule);
	}
}
