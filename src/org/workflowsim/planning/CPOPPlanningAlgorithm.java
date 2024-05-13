/**
 * Copyright 2024, Indian Institute of Technology (BHU), Varanasi 
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not
 * use this file except in compliance with the License. You may obtain a copy of
 * the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 * License for the specific language governing permissions and limitations under
 * the License.
 */
package org.workflowsim.planning;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.cloudbus.cloudsim.Consts;
import org.cloudbus.cloudsim.Log;
import org.workflowsim.CondorVM;
import org.workflowsim.FileItem;
import org.workflowsim.Task;
import org.workflowsim.utils.Parameters;
import java.util.PriorityQueue;

/**
 * The CPOP planning algorithm.
 * @author Bhanu Verma and Devashish Raj
 * @date April 15, 2024
 */


public class CPOPPlanningAlgorithm extends BasePlanningAlgorithm {
	/**
	 * These are the necessary variables and attributes to implement the algorithm.
	 * Added upwardRank and downwardRank in addition to the attributes of HEFT.
	 */
    private Map<Task, Map<CondorVM, Double>> computationCosts;	
    private Map<Task, Map<Task, Double>> transferCosts;		
    private Map<Task, Double> rank;	
    private Map<Task, Double> upwardRank;
    private Map<Task, Double> downwardRank;
    private Map<CondorVM, List<Event>> schedules;
    private Map<Task, Double> earliestFinishTimes;
    private double averageBandwidth;
    
    
    /**
     * trying to add a map to contain average computation cost for each task as described in algorithm
     * trying to get the entry task
     */
    private Map<Task, Double> avgComputaionCosts;
    private Task eTask;
    private List<Task> criticalPath;
    private CondorVM CP_VM;
    
    private class Event {

        public double start;
        public double finish;

        public Event(double start, double finish) {
            this.start = start;
            this.finish = finish;
        }
    }

    private class TaskRank implements Comparable<TaskRank> {

        public Task task;
        public Double rank;

        public TaskRank(Task task, Double rank) {
            this.task = task;
            this.rank = rank;
        }

        @Override
        public int compareTo(TaskRank o) {
            return o.rank.compareTo(rank);
        }
    }
    
    /**
     * Initialized all the variables to empty lists.
     */
    public CPOPPlanningAlgorithm() {
        computationCosts = new HashMap<>();
        transferCosts = new HashMap<>();
        rank = new HashMap<>();
        upwardRank = new HashMap<>();
        downwardRank = new HashMap<>();
        earliestFinishTimes = new HashMap<>();
        schedules = new HashMap<>();
        avgComputaionCosts = new HashMap<>();
        criticalPath = new ArrayList<Task>();
        CP_VM = null;
        eTask = null;
    }

    /**
     * The main function
     */
    @Override
    public void run() {
        Log.printLine("CPOP planner running with " + getTaskList().size()
                + " tasks and " + getVmList().size() + "vms.");
        System.out.println("Flow reaches here atleast");
        averageBandwidth = calculateAverageBandwidth();

        /**
         * saare Vms ke schedules ko empty list se initialize kr do
         */
        for (Object vmObject : getVmList()) {
            CondorVM vm = (CondorVM) vmObject;
            schedules.put(vm, new ArrayList<>());
        }

        // Prioritization phase
        calculateComputationCosts();
        calculateTransferCosts();
        calculateAvgComputationCosts();
        System.out.println("Why flow is here");
        calculateRanks();
        calculateEntryTask();
        calculateCriticalPath(eTask);
        
        /**
         * Printing for the sake for clarity --
         */
        System.out.println("Here entry task is " + eTask.getCloudletId());
        System.out.println("");
        System.out.println("Critical Path is:");
        for(Task t:criticalPath) {
        	System.out.print(t.getCloudletId() + "->");
        }
        System.out.println("");
        findCriticalProcessor();
        System.out.println("");
        
        // Selection phase
        allocateTasks();
    }

    /**
     * Calculates the average available bandwidth among all VMs in Mbit/s
     *
     * @return Average available bandwidth in Mbit/s
     * calculates as sum of bw of all the vms divided by number of Vms
     */
    private double calculateAverageBandwidth() {
    	System.out.println("here I came to calculate bandwidth");
        double avg = 0.0;
        for (Object vmObject : getVmList()) {
            CondorVM vm = (CondorVM) vmObject;
            avg += vm.getBw();
        }
        return avg / getVmList().size();
    }

    /**
     * trying to get the entry task.
     */
    private void calculateEntryTask() {
        for (Task task : getTaskList()) {
            if (task.getParentList() == null || task.getParentList().isEmpty()) {
                eTask = task;
                return;
            }
        }
        eTask = getTaskList().get(0);
        return;
    }

    private void calculateCriticalPath(Task curr) {
    	criticalPath.add(curr);
    	double max = 0.0;
    	boolean foundNextTask = false;
    	Task nextTask = getTaskList().get(0);
    	for(Task ch: curr.getChildList()) {
    		foundNextTask = true;
    		double r = rank.get(ch);
    		if(r>max) {
    			nextTask = ch;
    			max = r;
    		}
    	}
    	if(foundNextTask) {
    		calculateCriticalPath(nextTask);
    	}else {
    		return;
    	}
    }
    
    /**
     * selects the most efficient process for critical path --
     */
    public  void  findCriticalProcessor() {
    	double minfinishTime= Double.MAX_VALUE;
    	CondorVM  critical_processor = null;
    	CondorVM vm = null;
    	 for (Object vmObject : getVmList()) {
             vm = (CondorVM) vmObject;
    		 double finishTime=0.0;
    		  
	    	 for (Task t:getTaskList()) {
	    		 double minReadyTime = 0.0;
	    		 if(criticalPath.contains(t)) {
	    			 finishTime += findFinishTime(t, vm, minReadyTime, false); 
    	         }
    	     }
	    	 System.out.println("The finish time on " + vm.getId() + " VM is "+ finishTime);

	    	 if (finishTime < minfinishTime) {
	    		 minfinishTime=finishTime;
	    		 critical_processor = vm; 
	    	 }
	    }
    	CP_VM = critical_processor;
    	System.out.println("the Critical processor is "+ critical_processor.getId());
    	return;
    }
    /**
     * Populates the computationCosts field with the time in seconds to compute
     * a task in a vm.
     * 
     * Now we have the computation cost to run each task on each the virtual machine.
     */
    private void calculateComputationCosts() {
    	System.out.println("Came here to calculate computation costs.");
        for (Task task : getTaskList()) {
        	
        	System.out.println("Task:" + task.getCloudletId() + "length:" + task.getCloudletTotalLength());
            Map<CondorVM, Double> costsVm = new HashMap<>();
            for (Object vmObject : getVmList()) {
                CondorVM vm = (CondorVM) vmObject;
                if (vm.getNumberOfPes() < task.getNumberOfPes()) {
                    costsVm.put(vm, Double.MAX_VALUE);
                } else {
                    costsVm.put(vm,
                            task.getCloudletTotalLength() / vm.getMips());
                }
            }
            computationCosts.put(task, costsVm);
        }
    }
    
    private void calculateAvgComputationCosts() {
    	System.out.println("Came here");
    	for(Task task : getTaskList()) {
	    	double averageComputationCost = 0.0;
	
	        for (Double cost : computationCosts.get(task).values()) {
	            averageComputationCost += cost;
	        }
	
	        averageComputationCost /= computationCosts.get(task).size();
	        avgComputaionCosts.put(task,averageComputationCost);
    	}
    }

    /**
     * Populates the transferCosts map with the time in seconds to transfer all
     * files from each parent to each child
     */
    private void calculateTransferCosts() {
        // Initializing the matrix
        for (Task task1 : getTaskList()) {
            Map<Task, Double> taskTransferCosts = new HashMap<>();
            for (Task task2 : getTaskList()) {
                taskTransferCosts.put(task2, 0.0);
            }
            transferCosts.put(task1, taskTransferCosts);
        }

        // Calculating the actual values
        for (Task parent : getTaskList()) {
            for (Task child : parent.getChildList()) {
                transferCosts.get(parent).put(child,
                        calculateTransferCost(parent, child));
            }
        }
    }

    /**
     * Accounts the time in seconds necessary to transfer all files described
     * between parent and child
     *
     * @param parent
     * @param child
     * @return Transfer cost in seconds
     */
    private double calculateTransferCost(Task parent, Task child) {
        List<FileItem> parentFiles = parent.getFileList();
        List<FileItem> childFiles = child.getFileList();

        double acc = 0.0;

        for (FileItem parentFile : parentFiles) {
            if (parentFile.getType() != Parameters.FileType.OUTPUT) {
                continue;
            }

            for (FileItem childFile : childFiles) {
                if (childFile.getType() == Parameters.FileType.INPUT
                        && childFile.getName().equals(parentFile.getName())) {
                    acc += childFile.getSize();
                    break;
                }
            }
        }

        //file Size is in Bytes, acc in MB
        acc = acc / Consts.MILLION;
        // acc in MB, averageBandwidth in Mb/s
        return acc * 8 / averageBandwidth;
    }

    /**
     * Invokes calculateRank for each task to be scheduled
     */
    
    private void calculateRanks() {
    	System.out.println("I am being called but I am not being printed");
    	System.out.println();
    	System.out.println();
    	System.out.println("Task  	: 	Urank :		Drank :		Rank");
        for (Task task : getTaskList()) {
            double Ur = calculateUpwardRank(task);
            double Dr = calculateDownwardRank(task);
            double Sum = Ur+Dr;
            rank.put(task, Sum);
            System.out.println(task.getCloudletId() + " : " + Ur + " : " + Dr + " : " + Sum);
        }
        System.out.println();
    	System.out.println();
    }
   
    private double calculateUpwardRank(Task task) {
    	if (upwardRank.containsKey(task)) {
            return upwardRank.get(task);
        }

    	/**
    	 * double averageComputationCost = 0.0;
		*
        * for (Double cost : computationCosts.get(task).values()) {
            averageComputationCost += cost;
        * }

        * averageComputationCost /= computationCosts.get(task).size();
    	 */
        double averageComputationCost = avgComputaionCosts.get(task);

        double max = 0.0;
        for (Task child : task.getChildList()) {
            double childCost = transferCosts.get(task).get(child)
                    + calculateUpwardRank(child);
            max = Math.max(max, childCost);
        }

        upwardRank.put(task, averageComputationCost + max);

        return upwardRank.get(task);
    }
    
    private double calculateDownwardRank(Task task) {
    	if(downwardRank.containsKey(task)) {
    		return downwardRank.get(task);
    	}
    	double max = 0.0;
        for (Task parent : task.getParentList()) {
        	/**
        	 * double averageComputationCost = 0.0;

	            for (Double cost : computationCosts.get(parent).values()) {
	                averageComputationCost += cost;
	                
	            }

            	averageComputationCost /= computationCosts.get(task).size();
        	 */
        	double averageComputationCost = avgComputaionCosts.get(parent);
            if(task.getParentList() != null) {
            double parentCost = transferCosts.get(parent).get(task)
                    + calculateDownwardRank(parent)+ averageComputationCost ;
            max = Math.max(max, parentCost);}
            else max=0;
           
        }
        downwardRank.put(task, max);
        return downwardRank.get(task);
    }
    /**
     * Populates rank.get(task) with the rank of task as defined in the HEFT
     * paper.
     *
     * @param task The task have the rank calculates
     * @return The rank
     */
    
   /**
    * This code is not useful in CPOP--
    */
    /**
    private double calculateRank(Task task) {
        if (rank.containsKey(task)) {
            return rank.get(task);
        }

        double averageComputationCost = 0.0;

        for (Double cost : computationCosts.get(task).values()) {
            averageComputationCost += cost;
        }

        averageComputationCost /= computationCosts.get(task).size();

        double max = 0.0;
        for (Task child : task.getChildList()) {
            double childCost = transferCosts.get(task).get(child)
                    + calculateRank(child);
            max = Math.max(max, childCost);
        }

        rank.put(task, averageComputationCost + max);

        return rank.get(task);
    }
    */

    /**
     * Allocates all tasks to be scheduled in non-ascending order of schedule.
     */
    /**
    private void allocateTasks() {
    	for(Task t:getTaskList()) {
    		if(criticalPath.contains(t)) {
    			allocateCriticalTask(t);
    		}
    	}
        List<TaskRank> taskRank = new ArrayList<>();
        for (Task task : rank.keySet()) {
        	if(!criticalPath.contains(task)) {
        		taskRank.add(new TaskRank(task, rank.get(task)));
        	}
        }

        // Sorting in non-ascending order of rank
        Collections.sort(taskRank);
        for (TaskRank rank : taskRank) {
            allocateTask(rank.task);
        }

    }
    */
    
    private void allocateTasks() {
    	PriorityQueue<TaskRank> pq = new PriorityQueue<>();
    	for(Task task : rank.keySet()) {
    		if(task.getParentList()==null || task.getParentList().isEmpty()) {
    			pq.add(new TaskRank(task, rank.get(task)));
    		}
    	}
//    	pq.add(new TaskRank(eTask, rank.get(eTask)));
    	while (!pq.isEmpty()) {
            TaskRank t = pq.poll();
            Task currTask = t.task;
            System.out.println("Task ID: " + currTask.getCloudletId() + ", Rank: " + t.rank);
            if(criticalPath.contains(currTask)) {
            	allocateCriticalTask(currTask);
            	for(Task child : currTask.getChildList()) {
            		boolean isReady = true;
            		for(Task par: child.getParentList()) {
            			if(earliestFinishTimes.get(par)==null) {
            				isReady = false;
            				break;
            			}
            		}
            		if(isReady) {
            			pq.add(new TaskRank(child, rank.get(child)));
            		}
            	}
            }else {
            	allocateTask(currTask);
            	for(Task child : currTask.getChildList()) {
            		boolean isReady = true;
            		for(Task par: child.getParentList()) {
            			if(earliestFinishTimes.get(par)==null) {
            				isReady = false;
            				break;
            			}
            		}
            		if(isReady) {
            			pq.add(new TaskRank(child, rank.get(child)));
            		}
            	}
            }
        }
    }
    
    
    public void allocateCriticalTask(Task t) {
    	System.out.println("Allocating critical task " + t.getCloudletId() + " on critical VM " + CP_VM.getId());
    	CondorVM chosenVM = CP_VM;
    	double earliestFinishTime = Double.MAX_VALUE;
        double bestReadyTime = 0.0;
        double finishTime;
        double minReadyTime = 0.0;
        
        for (Task parent : t.getParentList()) {
            double readyTime = earliestFinishTimes.getOrDefault(parent, 0.0);
            if (parent.getVmId() != CP_VM.getId()) {
                readyTime += transferCosts.get(parent).get(t);
            }
            minReadyTime = Math.max(minReadyTime, readyTime);
        }
        finishTime = findFinishTime(t, CP_VM, minReadyTime, false);
        bestReadyTime = minReadyTime;
        earliestFinishTime = finishTime;
        chosenVM = CP_VM;
        findFinishTime(t, chosenVM, bestReadyTime, true);
        earliestFinishTimes.put(t, earliestFinishTime);
        t.setVmId(chosenVM.getId());
    }

    /**
     * Schedules the task given in one of the VMs minimizing the earliest finish
     * time
     *
     * @param task The task to be scheduled
     * @pre All parent tasks are already scheduled
     */
    private void allocateTask(Task task) {
    	System.out.print("Allocating task: " + task.getCloudletId());
        CondorVM chosenVM = null;
        double earliestFinishTime = Double.MAX_VALUE;
        double bestReadyTime = 0.0;
        double finishTime;

        for (Object vmObject : getVmList()) {
            CondorVM vm = (CondorVM) vmObject;
            double minReadyTime = 0.0;

            for (Task parent : task.getParentList()) {
                double readyTime = earliestFinishTimes.get(parent);
                if (parent.getVmId() != vm.getId()) {
                    readyTime += transferCosts.get(parent).get(task);
                }
                minReadyTime = Math.max(minReadyTime, readyTime);
            }

            finishTime = findFinishTime(task, vm, minReadyTime, false);

            if (finishTime < earliestFinishTime) {
                bestReadyTime = minReadyTime;
                earliestFinishTime = finishTime;
                chosenVM = vm;
            }
        }

        findFinishTime(task, chosenVM, bestReadyTime, true);
        earliestFinishTimes.put(task, earliestFinishTime);
        System.out.println(" On VM " + chosenVM.getId());
        task.setVmId(chosenVM.getId());
    }

    /**
     * Finds the best time slot available to minimize the finish time of the
     * given task in the vm with the constraint of not scheduling it before
     * readyTime. If occupySlot is true, reserves the time slot in the schedule.
     *
     * @param task The task to have the time slot reserved
     * @param vm The vm that will execute the task
     * @param readyTime The first moment that the task is available to be
     * scheduled
     * @param occupySlot If true, reserves the time slot in the schedule.
     * @return The minimal finish time of the task in the vmn
     */
    private double findFinishTime(Task task, CondorVM vm, double readyTime,
            boolean occupySlot) {
        List<Event> sched = schedules.get(vm);
        double computationCost = computationCosts.get(task).get(vm);
        double start, finish;
        int pos;

        if (sched.isEmpty()) {
            if (occupySlot) {
                sched.add(new Event(readyTime, readyTime + computationCost));
            }
            return readyTime + computationCost;
        }

        if (sched.size() == 1) {
            if (readyTime >= sched.get(0).finish) {
                pos = 1;
                start = readyTime;
            } else if (readyTime + computationCost <= sched.get(0).start) {
                pos = 0;
                start = readyTime;
            } else {
                pos = 1;
                start = sched.get(0).finish;
            }

            if (occupySlot) {
                sched.add(pos, new Event(start, start + computationCost));
            }
            return start + computationCost;
        }

        // Trivial case: Start after the latest task scheduled
        start = Math.max(readyTime, sched.get(sched.size() - 1).finish);
        finish = start + computationCost;
        int i = sched.size() - 1;
        int j = sched.size() - 2;
        pos = i + 1;
        while (j >= 0) {
            Event current = sched.get(i);
            Event previous = sched.get(j);

            if (readyTime > previous.finish) {
                if (readyTime + computationCost <= current.start) {
                    start = readyTime;
                    finish = readyTime + computationCost;
                }

                break;
            }
            if (previous.finish + computationCost <= current.start) {
                start = previous.finish;
                finish = previous.finish + computationCost;
                pos = i;
            }
            i--;
            j--;
        }

        if (readyTime + computationCost <= sched.get(0).start) {
            pos = 0;
            start = readyTime;

            if (occupySlot) {
                sched.add(pos, new Event(start, start + computationCost));
            }
            return start + computationCost;
        }
        if (occupySlot) {
            sched.add(pos, new Event(start, finish));
        }
        return finish;
    }
}


