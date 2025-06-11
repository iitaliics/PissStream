import pygame

"""
Basics:
 - Integration (apply external forces)
 - Projection (make fluid incompressable)
 - Advection (velocities move in the direction of velocity as if they were particles)

 https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/17-fluidSim.html#L94
"""

gravity = -9.81
over_relaxation = 1.9 # between 1 and 2
cnt = 0

FIELD = {
	'U_FIELD': 0,
	'V_FIELD': 1,
	'S_FIELD': 2,
}

class fluid:
    def __init__(self, density, xNum, yNum, cell_height):
        self.density = density
        self.xNum = xNum
        self.yNum = yNum
        self.num_cells = xNum * yNum
        self.height = cell_height

        self.u = [_ for _ in range(self.num_cells)]
        self.v = [_ for _ in range(self.num_cells)]

        self.new_u = [_ for _ in range(self.num_cells)]
        self.new_v = [_ for _ in range(self.num_cells)]

        self.p = [_ for _ in range(self.num_cells)]
        self.s = [_ for _ in range(self.num_cells)]  # Whether the cell is solid (0) or empty (1). Used to create barriers
        self.m = [1.0 for _ in range(self.num_cells)]
        self.new_m = [_ for _ in range(self.num_cells)]

    def integrate(self, dt):
            for x in range(1, self.xNum):
                for y in range(1, self.yNum - 1):
                    self.v[self.yNum * x + j] = gravity * dt

    def incompressibility_solver(self, iterations, dt):
        n = self.yNum #number of items per row for a flat array
        cp = self.density * self.height / dt #coefficient of pressure for simulation

        for _ in range(iterations):
            for i in range(1, self.xNum - 1):
                for j in range(1, self.yNum - 1):
                      
                    if self.s[n * i + j] == 0 : continue

                    sx0 = self.s[(i - 1)*n + j]
                    sx1 = self.s[(i + 1)*n + j]

                    sy0 = self.s[i*n + j - 1]
                    sy1 = self.s[i*n + j + 1]

                    s = sx0 + sx1 + sy0 + sy1
                    if s == 0 : continue #all surrounding sides are boundaries

                    divergence = self.u[(i + 1) * n + j] - self.u[i * n + j] + self.v[i * n + j + 1] - self.v[i * n + j]

                    pressure = -divergence / s

                    pressure *= over_relaxation
                    self.p[i * n + j] += cp * pressure

                    self.u[i * n + j] -= sx0 * pressure
                    self.u[(i + 1) * n + j] += sx1 * pressure
                    self.v[i * n + j] -= sy0 * pressure
                    self.v[i * n + j + 1] += sy1 * pressure

    def boundary_conditions(self):
        """
        Sets the edge cells to their closest inside neighbour
        """
        n = self.yNum #number of items per row

        for i in range(self.xNum):
            self.u[i * n + 0] = self.u[i * n + 1]
            self.u[i * n + self.yNum - 1] = self.u[i * n + self.yNum - 2]

        for j in range(self.yNum):
            self.v[0*n + j] = self.v[1 * n + j]
            self.v[(self.xNum - 1) * n + j] = self.v[(self.xNum - 2) * n + j] 

    def sample_field(self, x, y, field_type):
        n = self.yNum
        h = self.height
        h1 = 1 / h
        h2 = 0.5 * h

        x = max(min(x, self.xNum * h), h)
        y = max(min(y, self.xNum * h), h)

        match field_type:
            case 1: 
                f = self.u
                dy = h2

            case 2: 
                f = self.v 
                dy = h2

            case 3: 
                f = self.m 
                dx = h2
                dy = h2

        x0 = min(int((x - dx) * h1), self.xNum - 1)
        tx = ((x - dx) - x0 * h) * h1
        x1 = min(x0 + 1,  self.xNum - 1)

        y0 = min(int((y - dy) * h1), self.yNum - 1)
        ty = ((y - dy) - y0 * h) * h1
        y1 = min(y0 + 1,  self.yNum - 1)

        sx = 1.0 - tx
        sy = 1.0 - ty

        val = sx * sy * f[x0*n + y0] + \
            tx * sy * f[x1*n + y0] + \
            tx * ty * f[x1*n + y1] + \
            sx * ty * f[x0*n + y1]
        
        return val
    
    def avgU(self, i, j):
        n = self.yNum
        u = (self.u[i * n + j - 1] + self.u[i * n + j] + \
             self.u[(i + 1) * n + j - 1] + self.u[(i + 1) * n + j]) * 0.25
        return u
    
    def avgV(self, i, j):
        n = self.yNum
        v = (self.v[(i - 1) * n + j] + self.v[i * n + j] + \
             self.v[(i - 1) * n + j + 1] + self.v[i * n + j + 1]) * 0.25
        return v


    def advectVel(self, dt):
        self.newU = self.u
        self.newV = self.v
        n = self.yNum
        h = self.height
        h2 = 0.5 * h
        
        for i in range(1, self.xNum):
            for j in range(1, self.yNum):
    #			cnt++
                
                # u component
                if not self.s[i * n + j] == 0.0 and not self.s[(i - 1) * n + j] == 0.0 and j < self.yNum - 1:
                    x = i * h
                    y = j * h + h2
                    u = self.u[i * n + j]
                    v = self.avgV(i, j)
                    
                    x = x - dt * u
                    y = y - dt * v
                    u = self.sampleField(x, y, FIELD['U_FIELD'])
                    self.newU[i * n + j] = u
                    
                # v component
                if not self.s[i * n + j] == 0.0 and not self.s[(i - 1) * n + j - 1] == 0.0 and i < self.xNum - 1:
                    x = i * h + h2
                    y = j * h
                    u = self.avgU(i, j)
                    v = self.v[i * n + j]
                    
                    x = x - dt * u
                    y = y - dt * v
                    v = self.sampleField(x, y, FIELD['V_FIELD'])
                    self.newV[i * n + j] = v

    def advectSmoke(self, dt):
        self.newM = self.m
        
        n = self.yNum
        h = self.height
        h2 = 0.5 * h
        
        for i in range(1, self.xNum - 1):
            for j in range(1, self.yNum - 1):
                if not self.s[i * n + j] == 0.0:
                    u = (self.u[i * n + j] + self.u[(i + 1) * n + j]) * 0.5
                    v = (self.v[i * n + j] + self.v[i * n + j + 1]) * 0.5
                    x = i * h + h2 - dt * u
                    y = j * h + h2 - dt * v
                    
                    self.newM[i * n + j] = self.sampleField(x, y, FIELD['S_FIELD'])
        
        self.m = self.newM
        
    def simulate(self, dt, gravity, numIters):
        self.integrate(dt, gravity)
        
        self.p = [0.0 for _ in range(0, len(self.p))] 
        self.incompressibility_solver(numIters, dt)
        
        self.boundary_conditions()
        self.advectVel(dt)
        self.advectSmoke(dt)
      
tunnel = fluid(1000.0, 30, 30, 20)