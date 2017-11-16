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

/* ----------------------------------------------------------------------
   Contributing author: Mark Sears (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "plumedexchange.h"
#include "universe.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "integrate.h"
#include "modify.h"
#include "compute.h"
#include "force.h"
#include "output.h"
#include "thermo.h"
#include "fix.h"
#include "random_park.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define PLUMEDEXCHANGE_DEBUG 1

/* ---------------------------------------------------------------------- */

Plumed_Exchange::Plumed_Exchange(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

Plumed_Exchange::~Plumed_Exchange()
{
  MPI_Comm_free(&roots);
  if (ranswap) delete ranswap;
  delete ranboltz;
  delete [] set_ang;
  delete [] ang2world;
  delete [] world2ang;
  delete [] world2root;
}

/* ----------------------------------------------------------------------
   perform Plumed_Exchangeing with inter-world swaps
------------------------------------------------------------------------- */

void Plumed_Exchange::command(int narg, char **arg)
{
  if (universe->nworlds == 1)
    error->all(FLERR,"Must have more than one processor partition to Plumed_Exchange");
  if (domain->box_exist == 0)
    error->all(FLERR,"Plumed_Exchange command before simulation box is defined");
  if (narg != 8 && narg != 9)
    error->universe_all(FLERR,"Illegal Plumed_Exchange command");

  int nsteps = force->inumeric(FLERR,arg[0]);
  nevery = force->inumeric(FLERR,arg[1]);
  double ang = force->numeric(FLERR,arg[2]);

  for (whichfix = 0; whichfix < modify->nfix; whichfix++)
    if (strcmp(arg[3],modify->fix[whichfix]->id) == 0) break;
  if (whichfix == modify->nfix)
    error->universe_all(FLERR,"Plumed_Exchangeing fix ID is not defined");

  kappa = force->inumeric(FLERR,arg[4]);
  both_temp = force->inumeric(FLERR,arg[5]);
  seed_swap = force->inumeric(FLERR,arg[6]);
  seed_boltz = force->inumeric(FLERR,arg[7]);

  my_set_ang = universe->iworld;
  if (narg == 9) my_set_ang = force->inumeric(FLERR,arg[8]);

  // swap frequency must evenly divide total # of timesteps

  if (nevery == 0)
    error->universe_all(FLERR,"Invalid frequency in Plumed_Exchange command");
  nswaps = nsteps/nevery;
  if (nswaps*nevery != nsteps)
    error->universe_all(FLERR,"Non integer # of swaps in Plumed_Exchange command");

  // fix style must be appropriate for Plumed_Exchangeature control

//   printf("%s",modify->fix[whichfix]->style);
  if ((strcmp(modify->fix[whichfix]->style,"plumed") != 0))
    error->universe_all(FLERR,"Plumed_Exchangeing Plumed_Exchangeature fix is not valid");

//   fprintf("%f", modify->fix[whichfix]->kappa);

  // setup for long Plumed_Exchangeing run

  update->whichflag = 1;
  update->nsteps = nsteps;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + nsteps;
  if (update->laststep < 0 || update->laststep > MAXBIGINT)
    error->all(FLERR,"Too many timesteps");

  lmp->init();

  // local storage

  me_universe = universe->me;
  MPI_Comm_rank(world,&me);
  nworlds = universe->nworlds;
  iworld = universe->iworld;
  boltz = force->boltz;

  // pe_compute = ptr to thermo_pe compute
  // notify compute it will be called at first swap

  int id = modify->find_compute("thermo_pe");
  if (id < 0) error->all(FLERR,"Plumed_Exchangeing could not find thermo_pe compute");
  Compute *pe_compute = modify->compute[id];
  pe_compute->addstep(update->ntimestep + nevery);

  // create MPI communicator for root proc from each world

  int color;
  if (me == 0) color = 0;
  else color = 1;
  MPI_Comm_split(universe->uworld,color,0,&roots);

  // RNGs for swaps and Boltzmann test
  // warm up Boltzmann RNG

  if (seed_swap) ranswap = new RanPark(lmp,seed_swap);
  else ranswap = NULL;
  ranboltz = new RanPark(lmp,seed_boltz + me_universe);
  for (int i = 0; i < 100; i++) ranboltz->uniform();

  // world2root[i] = global proc that is root proc of world i

  world2root = new int[nworlds];
  if (me == 0)
    MPI_Allgather(&me_universe,1,MPI_INT,world2root,1,MPI_INT,roots);
  MPI_Bcast(world2root,nworlds,MPI_INT,0,world);

  // create static list of set Plumed_Exchange angles
  // allgather Plumed_Exchangeing arg "ang" across root procs
  // bcast from each root to other procs in world

  set_ang = new double[nworlds];
  if (me == 0) MPI_Allgather(&ang,1,MPI_DOUBLE,set_ang,1,MPI_DOUBLE,roots);
  MPI_Bcast(set_ang,nworlds,MPI_DOUBLE,0,world);

  // create world2ang only on root procs from my_set_ang
  // create ang2world on root procs from world2ang,
  //   then bcast to all procs within world

  world2ang = new int[nworlds];
  ang2world = new int[nworlds];
  if (me == 0) {
    MPI_Allgather(&my_set_ang,1,MPI_INT,world2ang,1,MPI_INT,roots);
    for (int i = 0; i < nworlds; i++) ang2world[world2ang[i]] = i;
  }
  MPI_Bcast(ang2world,nworlds,MPI_INT,0,world);

  // if restarting Plumed_Exchangeing, reset ang target of Fix to current my_set_ang

  if (narg == 7) {
    double new_ang = set_ang[my_set_ang];
    modify->fix[whichfix]->reset_target(new_ang);
  }

  // setup Plumed_Exchangeing runs

  int i,which,partner,swap,partner_set_ang,partner_world;
  double pe,pe_partner,boltz_factor,new_ang;
  MPI_Status status;
  // ***CHANGED***
  double curr_th, curr_th_partner = 0, centre, centre_partner;

  if (me_universe == 0 && universe->uscreen)
    fprintf(universe->uscreen,"Setting up Plumed_Exchangeing ...\n");

  update->integrate->setup();

  if (me_universe == 0) {
    if (universe->uscreen) {
      fprintf(universe->uscreen,"Step");
      for (int i = 0; i < nworlds; i++)
        fprintf(universe->uscreen," T%d",i);
      fprintf(universe->uscreen,"\n");
    }
    if (universe->ulogfile) {
      fprintf(universe->ulogfile,"Step");
      for (int i = 0; i < nworlds; i++)
        fprintf(universe->ulogfile," T%d",i);
      fprintf(universe->ulogfile,"\n");
    }
    print_status();
  }

  timer->init();
  timer->barrier_start(TIME_LOOP);

  for (int iswap = 0; iswap < nswaps; iswap++) {

    // run for nevery timesteps

    update->integrate->run(nevery);
//     double new_at = -2.0;
//     modify->fix[whichfix]->p->cmd("setAt", &new_at);
    modify->fix[whichfix]->reset_target(new_ang);

    // compute PE
    // notify compute it will be called at next swap

//     pe = pe_compute->compute_scalar();
//     pe_compute->addstep(update->ntimestep + nevery);
    // ***CHANGED***
    curr_th = modify->fix[whichfix]->compute_plumed_arg();
    // Checked this returns correct value

    // which = which of 2 kinds of swaps to do (0,1)

    if (!ranswap) which = iswap % 2;
    else if (ranswap->uniform() < 0.5) which = 0;
    else which = 1;

    // partner_set_ang = which set ang I am partnering with for this swap

    if (which == 0) {
      if (my_set_ang % 2 == 0) { partner_set_ang = my_set_ang + 1; }
      else { partner_set_ang = my_set_ang - 1; }
    } else {
      if (my_set_ang % 2 == 1) { partner_set_ang = my_set_ang + 1; }
      else { partner_set_ang = my_set_ang - 1; }
    }

    // partner = proc ID to swap with
    // if partner = -1, then I am not a proc that swaps

    partner = -1;
    if (me == 0 && partner_set_ang >= 0 && partner_set_ang < nworlds) {
      partner_world = ang2world[partner_set_ang];
      partner = world2root[partner_world];
    }

    // swap with a partner, only root procs in each world participate
    // hi proc sends PE to low proc
    // lo proc make Boltzmann decision on whether to swap
    // lo proc communicates decision back to hi proc

    // ***CHANGED***
    centre = set_ang[my_set_ang];
    centre_partner = set_ang[partner_set_ang];

    swap = 0;
    if (partner != -1) {
      // ***CHANGED***
      if (me_universe > partner) {
        MPI_Send(&curr_th,1,MPI_DOUBLE,partner,0,universe->uworld);
      }
      else {
        MPI_Recv(&curr_th_partner,1,MPI_DOUBLE,partner,0,universe->uworld,&status);
      }

//       printf("partner :%d, my curr angle: %f, partner curr angle: %f\n", partner, curr_th, curr_th_partner);

      // ***CHANGED***
      if (me_universe < partner) {
        double fac = (curr_th - curr_th_partner) * 
                     (centre - centre_partner) ;
        boltz_factor = 1.0 * kappa/(boltz*both_temp) * fac;
        if (boltz_factor >= 0.0) { swap = 1; }
        else if (ranboltz->uniform() < exp(boltz_factor)) { swap=1; }
      }

      if (me_universe < partner) {
        MPI_Send(&swap,1,MPI_INT,partner,0,universe->uworld); }
      else {
        MPI_Recv(&swap,1,MPI_INT,partner,0,universe->uworld,&status); }

#ifdef PLUMEDEXCHANGE_DEBUG
      if (me_universe < partner)
        printf("SWAP %d & %d: yes = %d,at angs = %g %g, curr angs = %g %g, Bz = %g %g\n",
               me_universe,partner,swap,centre,centre_partner,
               curr_th,curr_th_partner,boltz_factor,exp(boltz_factor));
#endif

    }

    // bcast swap result to other procs in my world

    MPI_Bcast(&swap,1,MPI_INT,0,world);

    // ***CHANGED***
    // velocity rescaling not necessary since we have the same temperature
//     // rescale kinetic energy via velocities if move is accepted
// 
//     if (swap) scale_velocities(partner_set_ang,my_set_ang);

    // if my world swapped, all procs in world reset ang target of Fix

    if (swap) {
      new_ang = set_ang[partner_set_ang];
      modify->fix[whichfix]->reset_target(new_ang);
    }

    // update my_set_ang and ang2world on every proc
    // root procs update their value if swap took place
    // allgather across root procs
    // bcast within my world

    if (swap) my_set_ang = partner_set_ang;
    if (me == 0) {
      MPI_Allgather(&my_set_ang,1,MPI_INT,world2ang,1,MPI_INT,roots);
      for (i = 0; i < nworlds; i++) ang2world[world2ang[i]] = i;
    }
    MPI_Bcast(ang2world,nworlds,MPI_INT,0,world);

    // print out current swap status

    if (me_universe == 0) print_status();
  }

  timer->barrier_stop(TIME_LOOP);

  update->integrate->cleanup();

  Finish finish(lmp);
  finish.end(1);

  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
}

/* ----------------------------------------------------------------------
   scale kinetic energy via velocities a la Sugita
------------------------------------------------------------------------- */

void Plumed_Exchange::scale_velocities(int t_partner, int t_me)
{
  double sfactor = sqrt(set_ang[t_partner]/set_ang[t_me]);

  double **v = atom->v;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    v[i][0] = v[i][0]*sfactor;
    v[i][1] = v[i][1]*sfactor;
    v[i][2] = v[i][2]*sfactor;
  }
}

/* ----------------------------------------------------------------------
   proc 0 prints current Plumed_Exchangeing status
------------------------------------------------------------------------- */

void Plumed_Exchange::print_status()
{
  if (universe->uscreen) {
    fprintf(universe->uscreen,BIGINT_FORMAT,update->ntimestep);
    for (int i = 0; i < nworlds; i++)
      fprintf(universe->uscreen," %d", world2ang[i]);
    fprintf(universe->uscreen,"\n");
  }
  if (universe->ulogfile) {
    fprintf(universe->ulogfile,BIGINT_FORMAT,update->ntimestep);
    for (int i = 0; i < nworlds; i++)
      fprintf(universe->ulogfile," %d",world2ang[i]);
    fprintf(universe->ulogfile,"\n");
    fflush(universe->ulogfile);
  }
}
