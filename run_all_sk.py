import minizinc
import time
import sys

def test_sk(start_rd, end_rd):
    model = minizinc.Model("sk_distinguisher.mzn")
    solver = minizinc.Solver.lookup("cp-sat")
    
    with open("results_sk_all.txt", "w") as f:
        f.write(f"Starting SK search for RD={start_rd} to {end_rd}\n")
        f.write("=========================================================\n\n")
    
    for rd in range(start_rd, end_rd + 1):
        print(f"[*] Testing SK RD = {rd} ...")
        sys.stdout.flush()
        inst = minizinc.Instance(solver, model)
        inst["RD"] = rd
        
        start_time = time.time()
        try:
            result = inst.solve()
            status = str(result.status)
            elapsed = time.time() - start_time
            
            with open("results_sk_all.txt", "a") as f:
                if "UNSATISFIABLE" in status:
                    f.write(f"[!] RD = {rd} is UNSATISFIABLE. (Took {elapsed:.2f}s)\n\n")
                elif "UNKNOWN" in status:
                    f.write(f"[!] RD = {rd} TIMED OUT. (Took {elapsed:.2f}s)\n\n")
                else:
                    f.write(f"[+] Found distinguisher for RD = {rd}! (Took {elapsed:.2f}s)\n")
                    f.write("="*80 + "\n")
                    f.write("FULL TRAIL DATA:\n")
                    f.write("="*80 + "\n")
                    f.write("--- UPPER TRAIL ---\n")
                    for r in range(rd + 1):
                        f.write(f"Round {r} AXU: {result['AXU'][r]}\n")
                        f.write(f"Round {r} DXU: {result['DXU'][r]}\n")
                    f.write("\n--- LOWER TRAIL ---\n")
                    for r in range(rd + 1):
                        f.write(f"Round {r} AXL: {result['AXL'][r]}\n")
                        f.write(f"Round {r} DXL: {result['DXL'][r]}\n")
                    f.write("="*80 + "\n\n")
        except Exception as e:
            with open("results_sk_all.txt", "a") as f:
                f.write(f"Error on SK RD={rd}: {e}\n\n")

def test_rk(setting_name, start_rd, end_rd, npt):
    model = minizinc.Model("distinguisherw_sat.mzn")
    solver = minizinc.Solver.lookup("cp-sat")
    
    with open("results_rk_all.txt", "a") as f:
        f.write(f"Starting {setting_name} (NPT={npt}) for RD={start_rd} to {end_rd}\n")
        f.write("=========================================================\n\n")
        
    for rd in range(start_rd, end_rd + 1):
        print(f"[*] Testing {setting_name} RD = {rd} ...")
        sys.stdout.flush()
        inst = minizinc.Instance(solver, model)
        inst["RD"] = rd
        inst["NPT"] = npt
        inst["is_related_tweakey"] = True
        
        start_time = time.time()
        try:
            result = inst.solve()
            status = str(result.status)
            elapsed = time.time() - start_time
            
            with open("results_rk_all.txt", "a") as f:
                if "UNSATISFIABLE" in status:
                    f.write(f"[!] {setting_name} RD = {rd} is UNSATISFIABLE. (Took {elapsed:.2f}s)\n\n")
                elif "UNKNOWN" in status:
                    f.write(f"[!] {setting_name} RD = {rd} TIMED OUT. (Took {elapsed:.2f}s)\n\n")
                else:
                    f.write(f"[+] Found {setting_name} distinguisher for RD = {rd}! (Took {elapsed:.2f}s)\n")
                    f.write("="*80 + "\n")
                    f.write("FULL TRAIL DATA:\n")
                    f.write("="*80 + "\n")
                    f.write(f"LANE (Related-Key Active Lanes): {result['LANE']}\n\n")
                    f.write("--- UPPER TRAIL ---\n")
                    for r in range(rd + 1):
                        f.write(f"Round {r} AXU: {result['AXU'][r]}\n")
                        f.write(f"Round {r} DXU: {result['DXU'][r]}\n")
                    f.write("\n--- LOWER TRAIL ---\n")
                    for r in range(rd + 1):
                        f.write(f"Round {r} AXL: {result['AXL'][r]}\n")
                        f.write(f"Round {r} DXL: {result['DXL'][r]}\n")
                    f.write("="*80 + "\n\n")
        except Exception as e:
            with open("results_rk_all.txt", "a") as f:
                f.write(f"Error on {setting_name} RD={rd}: {e}\n\n")

if __name__ == "__main__":
    print("Starting full overnight run...")
    
    # 1. Single-Key
    #test_sk(1, 10)
    
    # 2. Related-Key
    with open("results_rk_all.txt", "w") as f:
        f.write("uLBC Related-Key Impossible Differential Results\n\n")
        
    test_rk("Related-Key uLBC-128", 1, 14, 1)
    test_rk("Related-Key uLBC-256", 1, 14, 2)
    test_rk("Related-Key uLBC-384t", 1, 14, 3)
    
    print("All tasks completed.")
