package examples;

public class FitnessFunctionGWO {
	private static double [][] executiontimematrix, communicationtimematrix, taskoutputfilematrix, datatransfermatrix;
	private static double[] waittime;
	
	public int graph[][];
	int outputfilesize[];
	int mips[];
	double execcost[];
	double waitcost[];
	int tasklength[];
	int numOfTasks;
	int numOfVms;

	double commcost[][];
	FitnessFunctionGWO(double execcost[],double waitcost[],int mips[],int outputfilesize[],int tasklength[],int graph[][], int numOfTasks, int numOfVms)
	{
		this.tasklength = tasklength;
		this.execcost = execcost;
		this.mips = mips;
		this.outputfilesize = outputfilesize;
		this.graph = graph;
		this.waitcost = waitcost;
		this.numOfTasks = numOfTasks;
		this.numOfVms = numOfVms;
		commcost = new double[numOfTasks][numOfTasks];
		initializeMatrices();
	}
	
	private void initializeMatrices()
	{
		executiontimematrix = new double[numOfTasks][numOfVms];
		communicationtimematrix = new double[numOfTasks][numOfTasks];
		taskoutputfilematrix = new double[numOfTasks][numOfTasks];
		datatransfermatrix = new double[numOfTasks][numOfTasks];
		waittime = new double[numOfTasks];
			
		for(int i=0;i<numOfTasks;i++) {	
			for(int j=0;j<numOfVms;j++) {
				executiontimematrix[i][j] = tasklength[i]/mips[j];
				}
		}
		for(int i=0;i<numOfTasks;i++) {
			for(int j=i;j<numOfTasks;j++) {
				if(i==j)
					taskoutputfilematrix[i][j]=0;
				else
					taskoutputfilematrix[i][j]=outputfilesize[i]*graph[i][j];
			}
		}
		
		for(int i=0;i<numOfTasks;i++){
			for (int j=i;j<numOfTasks ;j++ ) {
				if(i==j)
					datatransfermatrix[i][j]=0;
				else{
					datatransfermatrix[i][j]=80;
					datatransfermatrix[j][i]=datatransfermatrix[i][j];
				}
			}
		}
		
		for(int i=0;i<numOfTasks;i++) {
			for(int j=i;j<numOfTasks;j++) {
				if(i==j)
					communicationtimematrix[i][j]=0;
				else
					{communicationtimematrix[i][j]=taskoutputfilematrix[i][j]/datatransfermatrix[i][j];
					}
				}
			}
		for(int i=0;i<numOfTasks;i++){
			for(int j=i;j<numOfTasks;j++){
				if(i==j)
					commcost[i][j]=0;
				else
				{
					commcost[i][j] = 3;
					commcost[j][i]=commcost[i][j];
				}
			}
		}
		printmatrices();
		
	}
	public double calculatecost(double[] position) {
		double cost = 0.0;
		double[] vmworkingcost = new double[numOfVms];
		
		for(int i=0;i<numOfTasks;i++){
			for(int j=i+1;j<numOfTasks;j++){
				if(taskoutputfilematrix[i][j]!=0) {
				waittime[j]=Math.max(waittime[j],waittime[i]+executiontimematrix[i][(int)position[i]]+communicationtimematrix[i][j]);
				}
			}
		}
		for(int i=0;i<numOfTasks;i++) {
			int vm = (int) position [i];
			vmworkingcost[vm]+=(executiontimematrix[i][vm])*execcost[vm];
		}
		for(int i=0;i<numOfTasks;i++) {
			int vm = (((int) position[i]));
			for(int j=i+1;j<numOfTasks;j++) {
			vmworkingcost[vm]+=(communicationtimematrix[i][j])*commcost[i][j];
			}
		}
		for(int i=0;i<numOfVms;i++)
			cost+=vmworkingcost[i];
		for(int i=0;i<numOfTasks;i++) {
			cost+= waittime[i]*waitcost[(int)position[i]];
		}
		return cost;
	}
	public double evaluate(double[] position) {
		return calculatecost(position);
	}
	public void printmatrices() {
		System.out.println("Execution Time Matrix");
		for(int i=0;i<numOfTasks;i++) {
			for(int j=0;j<numOfVms;j++) {
				System.out.print(executiontimematrix[i][j]+"\t");
			}
			System.out.println();
		}
		System.out.println("Communication Time Matrix");
		for(int i=0;i<numOfTasks;i++) {
			for(int j=0;j<numOfTasks;j++) {
				System.out.print(communicationtimematrix[i][j]+"\t");
			}
			System.out.println();
		}
		System.out.println("Communication cost Matrix");
		for(int i=0;i<numOfTasks;i++) {
			for(int j=0;j<numOfTasks;j++) {
				System.out.print(commcost[i][j]+"\t");
			}
			System.out.println();
		}
	}
	public double[][] getexecutiontimematrix(){
		return executiontimematrix;
	}
	public double[][] getcommunicationtimematrix(){
		return communicationtimematrix;
	}
	public double[][] getdatatransfermatrix(){
		return datatransfermatrix;
	}
	public double[][] getcommcost(){
		return commcost;
	}
	public double[][] gettaskoutputfilematrix(){
		return taskoutputfilematrix;
	}

}
