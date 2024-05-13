package org.workflowsim.scheduling;

import java.util.Iterator;
import org.cloudbus.cloudsim.Cloudlet;
import org.workflowsim.CondorVM;
import org.workflowsim.WorkflowSimTags;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import org.cloudbus.cloudsim.Vm;

public class SJFSchedulingAlgorithm extends BaseSchedulingAlgorithm {
	@Override
    public void run() {
        List<Cloudlet> cloudletList = getCloudletList();
        Collections.sort(cloudletList, new Comparator<Cloudlet>() {
            @Override
            public int compare(Cloudlet c1, Cloudlet c2) {
                return Long.compare(c1.getCloudletLength(), c2.getCloudletLength());
            }
        });
        for (Iterator<Cloudlet> it = cloudletList.iterator(); it.hasNext();) {
            Cloudlet cloudlet = it.next();
            System.out.println(cloudlet.getCloudletLength());
            boolean stillHasVm = false;
            for (Iterator<Vm> itc = getVmList().iterator(); itc.hasNext();) {

                CondorVM vm = (CondorVM) itc.next();
                if (vm.getState() == WorkflowSimTags.VM_STATUS_IDLE) {
                    stillHasVm = true;
                    vm.setState(WorkflowSimTags.VM_STATUS_BUSY);
                    cloudlet.setVmId(vm.getId());
                    System.out.println("vm"+vm.getId()+".mips: "+vm.getMips()+"  host: "+vm.getHost().getId());
                    getScheduledList().add(cloudlet);
                    break;
                }
            }
           
            if (!stillHasVm) {
                break;
            }
        }
    }
}
