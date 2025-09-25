#!/usr/bin/env python3
# Инженерийн гидравлик: Дөхүүлэх сувгийн гидравлик тооцоо (FORTRAN -> Python)
# Эх код (FORTRAN): Б.Аюурзана, 2017  |  Python хөрвүүлэлт (fixed): 2025

import math
import argparse
import sys
import shutil
import subprocess

def compute_channel(q, sl, mn, nb, h, jcz, b0, n=100000, alh=1e-4):
    """
    Core computation. Iterates bed width b from b0 with step alh until qt >= 1.2*q or n steps.
    Returns dict with results and writes 'qfb.dat' (columns: b, qt).
    """
    if nb <= 0.0:
        raise ValueError("Barzgarshilt (nb) must be > 0.")
    if h <= 0.0:
        raise ValueError("Depth h must be > 0.")
    if sl < 0.0:
        raise ValueError("Slope sl must be >= 0.")
    if jcz not in (0, 1, 2):
        raise ValueError("jcz must be 0 (Pawlowski), 1 (Agroskin), or 2 (Manning).")

    b = float(b0)
    bt = 0.0          # best/last width where qt <= q
    qth = 0.0
    rgt = 0.0
    czt = 0.0

    with open("qfb.dat", "w") as f:
        for _ in range(int(n)):
            # Geometry
            w = (b + mn * h) * h
            ksi = b + 2.0 * h * math.sqrt(1.0 + mn ** 2)
            if ksi <= 0.0:
                ksi = 1e-12
            rg = w / ksi  # hydraulic radius

            # Chezy C
            if jcz == 0:  # Pawlowski
                ay = 1.5 * math.sqrt(nb) if rg >= 1.0 else 1.3 * math.sqrt(nb)
                cz = (1.0 / nb) * (rg ** ay)
            elif jcz == 1:  # Agroskin
                # safe log10
                cz = (1.0 / nb) + 17.72 * math.log10(rg if rg > 1e-12 else 1e-12)
            else:  # Manning (Chezy via Manning)
                cz = (1.0 / nb) * (rg ** (1.0 / 6.0))

            # Discharge
            rgsl = rg * sl
            qt = w * cz * math.sqrt(rgsl if rgsl > 0.0 else 0.0)

            # Write pair (b, qt)
            f.write(f"{b:.10f} {qt:.10f}\n")

            # Track last value not exceeding target q
            if qt <= q:
                bt = b
                qth = qt
                rgt = rg
                czt = cz

            # Stop once we exceed 1.2*q
            if qt < 1.2 * q:
                b += alh
            else:
                break

    return {
        "bt": bt, "qth": qth, "rgt": rgt, "czt": czt,
        "cr05": czt * math.sqrt(rgt if rgt > 0.0 else 0.0)
    }

def parse_args():
    p = argparse.ArgumentParser(description="Offtake channel hydraulic calculation (Chezy options).")
    p.add_argument("--q", type=float, help="Discharge (m^3/s)")
    p.add_argument("--sl", type=float, help="Slope (m/m)")
    p.add_argument("--mn", type=float, help="Side slope (z : 1)")
    p.add_argument("--nb", type=float, help="Roughness parameter n or similar (dimensionless)")
    p.add_argument("--h", type=float, help="Depth h (m)")
    p.add_argument("--jcz", type=int, choices=[0, 1, 2],
                   help="Chezy method: 0=Pawlowski, 1=Agroskin, 2=Manning")
    p.add_argument("--b0", type=float, help="Initial bottom width b0 (m)")
    p.add_argument("--n", type=int, default=100000, help="Max iterations (default: 100000)")
    p.add_argument("--alh", type=float, default=1e-4, help="Step size for b (default: 1e-4 m)")
    p.add_argument("--plot", action="store_true", help="Run gnuplot with qfb.plt if available")
    p.add_argument("--noplot", action="store_true", help="Do not run gnuplot")
    return p.parse_args()

def prompt_float(msg):
    while True:
        s = input(msg).strip()
        try:
            return float(s)
        except ValueError:
            print("Please enter a number.")

def prompt_int(msg, allowed=None):
    while True:
        s = input(msg).strip()
        try:
            v = int(s)
            if allowed is None or v in allowed:
                return v
            print("Please enter one of:", allowed)
        except ValueError:
            print("Please enter an integer.")

def main():
    args = parse_args()

    # Interactive if any required argument is missing
    interactive = any(v is None for v in (args.q, args.sl, args.mn, args.nb, args.h, args.jcz, args.b0))

    if interactive:
        print("Zartsuulga ugnu uu, m3s-1")
        q = prompt_float("> ")
        print("Hewgii, Hajuu naluu, Barzgariig ugnu uu (3 too, зайгаар тусгаарлаад оруулна)")
        sl = prompt_float("  Hewgii (sl): ")
        mn = prompt_float("  Hajuu naluu (mn): ")
        nb = prompt_float("  Barzgarshilt (nb): ")
        print("h-g ugnu uu, m")
        h = prompt_float("> ")
        print("Shezi bodson arga? 0-Pawlowski, 1-Agroskin, 2-Manning")
        jcz = prompt_int("> ", allowed={0, 1, 2})
        print("Ehnii b ugnu uu, m")
        b0 = prompt_float("> ")
        n = args.n
        alh = args.alh
        do_plot = args.plot and not args.noplot
    else:
        q, sl, mn, nb, h, jcz, b0 = args.q, args.sl, args.mn, args.nb, args.h, args.jcz, args.b0
        n = args.n
        alh = args.alh
        do_plot = args.plot and not args.noplot

    try:
        res = compute_channel(q, sl, mn, nb, h, jcz, b0, n=n, alh=alh)
    except Exception as e:
        print("Error:", e)
        sys.exit(1)

    print("Tootsoot urgun ba zartsuulga: {:.6f} m, {:.6f} m3/s".format(res["bt"], res["qth"]))
    print("Tootsoot R ba C, CR^0.5: R={:.6f} m, C={:.6f}, C*sqrt(R)={:.6f}".format(
        res["rgt"], res["czt"], res["cr05"]
    ))
    print("Data written to 'qfb.dat' (columns: b, qt).")

    # Optional: call gnuplot if requested and available
    if do_plot:
        if shutil.which("gnuplot") is not None:
            try:
                subprocess.run(["gnuplot", "-persist", "qfb.plt"], check=False)
            except Exception as e:
                print("Note: gnuplot run error:", e)
        else:
            print("Note: gnuplot not found. Install it or run without --plot.")

if __name__ == "__main__":
    main()
