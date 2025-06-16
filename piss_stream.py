import pygame

"""
Basics:
 - Integration (apply external forces)
 - Projection (make fluid incompressable)
 - Advection (velocities move in the direction of velocity as if they were particles)

 https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/17-fluidSim.html#L94
"""

gravity = -9.81
over_relaxation = 1.9 # between 1 and 2 1.9
cnt = 0

FIELD = {
	'U_FIELD': 0,
	'V_FIELD': 1,
	'S_FIELD': 2,
}

class fluid:
    def __init__(self, density, xNum, yNum, cell_height):
        self.density = density
        self.xNum = xNum + 2
        self.yNum = yNum + 2
        self.num_cells = self.xNum * self.yNum 
        self.height = cell_height

        self.u = [0 for _ in range(self.num_cells)]
        self.v = [0 for _ in range(self.num_cells)]

        self.new_u = [0 for _ in range(self.num_cells)]
        self.new_v = [0 for _ in range(self.num_cells)]

        self.p = [0 for _ in range(self.num_cells)]
        self.s = [1 for _ in range(self.num_cells)]  # Whether the cell is solid (0) or empty (1). Used to create barriers
        self.m = [1.0 for _ in range(self.num_cells)]
        self.new_m = [0 for _ in range(self.num_cells)]

    def integrate(self, dt, grav):
            for x in range(1, self.xNum):
                for y in range(1, self.yNum - 1):
                    if self.s[x * self.yNum + y] == 0 and self.s[x * self.yNum + y - 1 == 0]:
                        self.v[self.yNum * x + y] += grav * dt

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

                    # if i == 15 and j == 15:
                    #     print("divergence =", divergence)

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

        # for i in range(1, self.xNum - 1):
        #     for j in range(1, self.yNum - 1):
        #         num = i * self.yNum + j
        #         if self.s[num] == 0:
        #             self.u[num] = 0
        #             self.v[num] = 0
        #             self.u[(i + 1) * self.yNum + j] = 0
        #             self.v[i * self.yNum + j + 1] = 0

    def sample_field(self, x, y, field_type):
        n = self.yNum
        h = self.height
        h1 = 1.0 / h
        h2 = 0.5 * h

        x = max(min(x, self.xNum * h), h)
        y = max(min(y, self.yNum * h), h)

        dx, dy = 0.0, 0.0
        f = []

        match field_type:
            case 0: 
                f = self.u.copy()
                dx = h2

            case 1: 
                f = self.v.copy()
                dy = h2

            case 2: 
                f = self.m.copy()
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
        self.newU = self.u.copy()
        self.newV = self.v.copy()
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
                    u = self.sample_field(x, y, FIELD['U_FIELD'])
                    self.newU[i * n + j] = u
                    
                # v component
                if not self.s[i * n + j] == 0.0 and not self.s[i * n + j - 1] == 0.0 and i < self.xNum - 1:
                    x = i * h + h2
                    y = j * h
                    u = self.avgU(i, j)
                    v = self.v[i * n + j]
                    
                    x = x - dt * u
                    y = y - dt * v
                    v = self.sample_field(x, y, FIELD['V_FIELD'])
                    self.newV[i * n + j] = v

        self.u = self.newU.copy()
        self.v = self.newV.copy()


    def advectSmoke(self, dt):
        self.newM = self.m.copy()
        
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
                    
                    self.newM[i * n + j] = self.sample_field(x, y, FIELD['S_FIELD'])
            
        self.m = self.newM.copy()

    def setup_walls(self):
        for y in range(self.yNum):      # ver
            for x in range(self.xNum):  # hor
                num = x * self.yNum + y 
                if y == 0 or x == 0 or y == self.yNum - 1:
                    self.s[num] = 0.0
                else:
                    self.s[num] = 1.0

    def create_obstable(self, xPos, yPos, radius):
        for x in range(self.xNum):
            for y in range(self.yNum):
                if (pow(xPos - x, 2) + pow(yPos - y, 2) < pow(radius, 2)):
                    self.s[(x) * self.yNum + y] = 0
                    self.s[(x) * self.yNum + y] = 0
            


    def wind_tunnel_behaviour(self, vel):
        # Add the wind velocity

        
        # for i in range(self.xNum - 1):
        #     for j in range(self.yNum - 1):
        #         if self.s[i * self.yNum + j] == 0:
        #             self.u[i * self.yNum + j] = 0.0
        #             self.u[(i + 1) * self.yNum + j] = 0.0
        #             self.v[i * self.yNum + j] = 0.0
        #             self.v[i * self.yNum + j + 1] = 0.0

        for j in range (1, self.yNum - 1):
            self.u[1 * self.yNum + j] = vel

        # Add the smoke ""mass"" to make it visible and see how it distributes
        
        for j in range (int(2 * self.yNum / 5), int(3 * self.yNum / 5)):
            self.m[j] = 21.0


    
        
    def simulate(self, dt, gravity, numIters, windSpeed):

        self.integrate(dt, gravity)
        
        self.p = [0.0 for _ in range(0, self.num_cells)] 
        self.incompressibility_solver(numIters, dt)
        
        self.boundary_conditions()
        self.advectVel(dt)
        self.advectSmoke(dt)

        self.wind_tunnel_behaviour(windSpeed)
        
def get_sci_color(val, min_val, max_val):
    # Clamp value between min_val and max_val - 0.0001
    val = min(max(val, min_val), max_val - 0.0001)

    d = max_val - min_val
    val = 0.5 if d == 0.0 else (val - min_val) / d

    m = 0.25
    num = int(val / m)
    s = (val - num * m) / m

    r, g, b = 0.0, 0.0, 0.0

    if num == 0:
        r, g, b = 0.0, s, 1.0
    elif num == 1:
        r, g, b = 0.0, 1.0, 1.0 - s
    elif num == 2:
        r, g, b = s, 1.0, 0.0
    elif num == 3:
        r, g, b = 1.0, 1.0 - s, 0.0

    return (int(255 * r), int(255 * g), int(255 * b))


      

window_width = 50
window_height = 25

cell_height = 3

tunnel = fluid(10.0, window_width, window_height, cell_height)
tunnel.setup_walls()
tunnel.create_obstable(15, int(window_height / 2) + 2, 4)

pygame.init()
screen = pygame.display.set_mode(((window_width * cell_height), (window_height * cell_height)))

done = False
while not (done):
    for event in pygame.event.get():  # User did something
        if event.type == pygame.QUIT:  # If user clicked close
            done = True  # Flag that we are done so we exit this loop


    
    gravity = 9.81
    gravity = 0
    # for _ in range(12):
        # tunnel.p[_] = 0.0
    tunnel.simulate(0.001, gravity, 100, 500)
    # tunnel.simulate(1.0 / 120, 9.81, 100, )

    

    

    # for j in range (tunnel.yNum):
    #     tunnel.m[tunnel.yNum + 1] = 1.0

    # for i in range(tunnel.xNum):
    #     tunnel.m[i * tunnel.yNum + 4] = 1.0
        
    
    # print(tunnel.s)
    # print(tunnel.p)
    # print(tunnel.m)
    # print(tunnel.s)

    minP = 0
    maxP = 0

    minM = 0
    maxM = 0
    
    for i in range(tunnel.num_cells):
        minP = min(minP, tunnel.p[i])
        maxP = max(maxP, tunnel.p[i])

    maxP = maxP if maxP > 0 else 0.000000000000001

    for i in range(tunnel.num_cells):
        minM = min(minM, tunnel.m[i])
        maxM = max(maxM, tunnel.m[i]) + 0.0005

    # print(maxM, minM)
    print(maxP, minP)
    for _ in range(tunnel.xNum):
        print(_, tunnel.p[_ * tunnel.yNum + 15])
        print(_, tunnel.p[_ * tunnel.yNum + 5])
    screen.fill((0, 0, 0))
    
    # debug
    # item = 15
    # print(tunnel.s[(item + 1) * tunnel.yNum + 15])
    # print(tunnel.m[(item + 1) * tunnel.yNum + 15])
    # print(tunnel.p[(item + 1) * tunnel.yNum + 15])
    # print(tunnel.u[(item + 1) * tunnel.yNum + 15])
    # print(tunnel.v[(item + 1) * tunnel.yNum + 15])


    for y in range(1, tunnel.yNum - 1):      # ver
        for x in range(1, tunnel.xNum - 1):  # hor
            num = x * tunnel.yNum + y  # match fluid class

            valM = 255 * ((tunnel.m[num] - minM) / (maxM - minM))
            # valP = 255 * ((tunnel.m[num] - minP) / (maxP - minP))

            colour = get_sci_color(tunnel.p[num], minP, maxP)


            # print(val)

            # print(valP, maxP, tunnel.m[num])

            # r, g, b = 0, valM, valP
            # r, g, b, = valM, valM, valM

            r = colour[0] - valM
            g = colour[1] - valM
            b = colour[2] - valM

            # r = colour[0] - (255 - valM)
            # g = colour[1] - (255 - valM)
            # b = colour[2] - (255 - valM)

            # r = colour[0]
            # g = colour[1] 
            # b = colour[2] 

            r = max(0, min(r, 255))
            b = max(0, min(b, 255))
            g = max(0, min(g, 255))

            if tunnel.s[num] == 0:
                r, g, b = 255, 0, 0

            # If you want to skip the 1-cell buffer borders:
            if 1 <= x < tunnel.xNum - 1 and 1 <= y < tunnel.yNum - 1:
                pygame.draw.rect(
                    screen,
                    (r, g, b),
                    [(x - 1) * cell_height, (y - 1) * cell_height, cell_height, cell_height]
                )

    pygame.display.update()

    