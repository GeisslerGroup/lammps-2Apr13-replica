/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(plumed_exchange,Plumed_Exchange)

#else

#ifndef LMP_PLUMED_EXCHANGE_H
#define LMP_PLUMED_EXCHANGE_H

#include "pointers.h"

namespace LAMMPS_NS {

class Plumed_Exchange : protected Pointers {
 public:
  Plumed_Exchange(class LAMMPS *);
  ~Plumed_Exchange();
  void command(int, char **);

 private:
  int me,me_universe;          // my proc ID in world and universe
  int iworld,nworlds;          // world info
  double boltz;                // copy from output->boltz
  MPI_Comm roots;              // MPI comm with 1 root proc from each world
  class RanPark *ranswap,*ranboltz;  // RNGs for swapping and Boltz factor
  int nevery;                  // # of timesteps between swaps
  int nswaps;                  // # of Plumed_Exchangeing swaps to perform
  double kappa;		       // spring constant for simulations
  int seed_swap;               // 0 = toggle swaps, n = RNG for swap direction
  int seed_boltz;              // seed for Boltz factor comparison
  int whichfix;                // index of Plumed_Exchangeature fix to use
  int fixstyle;                // what kind of Plumed_Exchangeature fix is used

  int my_set_ang;             // which set ang I am simulating
  double *set_ang;            // static list of replica set Plumed_Exchangeatures
  int *ang2world;             // ang2world[i] = world simulating set ang i
  int *world2ang;             // world2ang[i] = ang simulated by world i
  int *world2root;             // world2root[i] = root proc of world i

  void scale_velocities(int, int);
  void print_status();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Must have more than one processor partition to Plumed_Exchange

Cannot use the Plumed_Exchange command with only one processor partition.  Use
the -partition command-line option.

E: Plumed_Exchange command before simulation box is defined

The Plumed_Exchange command cannot be used before a read_data, read_restart, or
create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Plumed_Exchangeing fix ID is not defined

The fix ID specified by the Plumed_Exchange command does not exist.

E: Invalid frequency in Plumed_Exchange command

Nevery must be > 0.

E: Non integer # of swaps in Plumed_Exchange command

Swap frequency in Plumed_Exchange command must evenly divide the total # of
timesteps.

E: Plumed_Exchangeing Plumed_Exchangeature fix is not valid

The fix specified by the Plumed_Exchange command is not one that controls
ang (nvt or langevin).

E: Too many timesteps

The cummulative timesteps must fit in a 64-bit integer.

E: Plumed_Exchangeing could not find thermo_pe compute

This compute is created by the thermo command.  It must have been
explicitly deleted by a uncompute command.

*/
