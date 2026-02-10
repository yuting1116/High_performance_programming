#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "graphics.h"

// graphic +++
const float circleRadius=0.005, circleColor=0;
const int windowWidth=800;

void keep_within_box(float* xA, float* yA) {
  if(*xA > 1)
    *xA = 0;
  if(*yA > 1)
    *yA = 0;
}

// graphic ---

typedef struct {
    double pos_x;
    double pos_y;
    double mass;
    double vel_x;
    double vel_y;
    double brightness;
} Particle;

Particle* read_configuration(char* fileName, int N)
{
    FILE *fp;
    fp = fopen(fileName, "r");

    if (fp == NULL) {
        printf("open file error\n");
        return NULL;
    }

    Particle* particles = (Particle*)malloc(N * sizeof(Particle));

    for (int i = 0; i < N; i++)
    {   
        if (fread(&particles[i], sizeof(Particle), 1, fp) != 1)
        {
            printf("fread failed at particle %d\n", i);
            exit(EXIT_FAILURE);
        }
        // printf("Read file: Particel%d: mass(%.2f) pos(%.2f, %.2f) vel(%.2f, %.2f) brightness(%.2f)\n", i, particles[i].mass,
        //     particles[i].pos_x, particles[i].pos_y, particles[i].vel_x, particles[i].vel_y, particles[i].brightness);
    }

    fclose(fp);
    return particles;
}

void update_particles(int N, Particle*  particles, Particle*  particles_next, double delta_t)
{
    double G = 100.0f/N;
    double epsilon = 1e-3;
    double* f_x_ij = (double*)calloc(N, sizeof(double));
    double* f_y_ij = (double*)calloc(N, sizeof(double));
    for (int i = 0; i < N; i++)
    {        
        //Copy particle[i] to local variables avoiding reading memory each time
        double pos_x_i = particles[i].pos_x; 
        double pos_y_i = particles[i].pos_y;
        double vel_x_i = particles[i].vel_x;
        double vel_y_i = particles[i].vel_y;
        double mass_i = particles[i].mass;

        // j starts from i+1 to calculate each particle pair only once
        for (int j = i+1; j < N; j++)
        {
            double mass_j = particles[j].mass;
            double pos_x_j = particles[j].pos_x;
            double pos_y_j = particles[j].pos_y;
    
            double r_x = pos_x_i - pos_x_j;
            double r_y = pos_y_i - pos_y_j;
            double r_ij_mag = sqrt(r_x*r_x + r_y*r_y);
            double div = 1/((r_ij_mag + epsilon)*(r_ij_mag + epsilon)*(r_ij_mag + epsilon));
            // Update forces symmetrically for both particles: F_ij = -F_ji
            f_x_ij[i] += mass_j * div * r_x;
            f_y_ij[i] += mass_j * div * r_y;
            f_x_ij[j] += - mass_i * div * r_x;
            f_y_ij[j] += - mass_i * div * r_y;
        }

        double F_x = -G * mass_i * f_x_ij[i];
        double F_y = -G * mass_i * f_y_ij[i];
        double mass_i_div = 1.0f / mass_i;
        double vel_x_next = vel_x_i + delta_t * F_x * mass_i_div;
        double vel_y_next = vel_y_i + delta_t * F_y * mass_i_div;
        double pos_x_next = pos_x_i + delta_t * vel_x_next;
        double pos_y_next = pos_y_i + delta_t * vel_y_next;

        // update
        particles_next[i].pos_x = pos_x_next;
        particles_next[i].pos_y = pos_y_next;
        particles_next[i].vel_x = vel_x_next;
        particles_next[i].vel_y = vel_y_next;
        // printf("Particel%d: mass(%.2f) pos(%.2f, %.2f) vel(%.2f, %.2f) brightness(%.2f)\n", i, mass_i,
        //     particles_next[i].pos_x, particles_next[i].pos_y, particles_next[i].vel_x, particles_next[i].vel_y, particles_next[i].brightness);
    }

    free(f_x_ij);
    free(f_y_ij);
    f_x_ij = NULL;
    f_y_ij = NULL;
}

void write_history(FILE* fp, int N, Particle* particles)
{
    for (int i = 0; i < N; i++)
    {
        fwrite(&particles[i].pos_x, sizeof(double), 1, fp);
        fwrite(&particles[i].pos_y, sizeof(double), 1, fp);
        fwrite(&particles[i].mass, sizeof(double), 1, fp);
        fwrite(&particles[i].vel_x, sizeof(double), 1, fp);
        fwrite(&particles[i].vel_y, sizeof(double), 1, fp);
        fwrite(&particles[i].brightness, sizeof(double), 1, fp);
    }
    fclose(fp);
}

int main(int argc, char *argv[]) {
    
    if (argc <= 5)
    {
        printf("Not enough input!\n");
        return 0;
    }

    int N  = atoi(argv[1]);
    char* filename = argv[2];
    int nsteps = atoi(argv[3]);
    double delta_t = strtod(argv[4], NULL);
    bool en_graphics = (strcmp(argv[5], "1") == 0);
    printf("N:%d, filename:%s nstep:%d delta_t:%.2f, graphic:%d\n", N, filename, nsteps, delta_t, en_graphics);

    FILE *fp = fopen("result.gal", "wb");

    Particle* particles = read_configuration(filename, N);
    Particle* particles_next = (Particle*)malloc(N * sizeof(Particle));
    memcpy(particles_next, particles, N * sizeof(Particle));

    if (en_graphics)
    {
        float L=1, W=1;
        InitializeGraphics(argv[0], windowWidth, windowWidth/2);
        SetCAxes(0,1);

        for (int step = 0; step < nsteps; step++)
        {
            ClearScreen();
            for (int i = 0; i < N; i++)
            {
                DrawCircle(particles[i].pos_x, particles[i].pos_y, L, W, circleRadius, 1.0/(N+2)*(i+1));
            }

            update_particles(N, particles, particles_next, delta_t);
            memcpy(particles, particles_next, N * sizeof(Particle));

            Refresh();
            usleep(10000);
        }
        FlushDisplay();
        CloseDisplay();
    }
    else
    {
        for (int step = 0; step < nsteps; step++)
        {
            update_particles(N, particles, particles_next, delta_t);
            memcpy(particles, particles_next, N * sizeof(Particle));
        }
    }
    
    write_history(fp, N, particles);
    
    // To do
    // free struct Particle
    free(particles);
    free(particles_next);

    particles = NULL;
    particles_next = NULL;
    return 0;

}
