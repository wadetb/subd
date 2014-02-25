SubD Notes
==========

# TODO

- Don't build a QE mesh for the final level.
- Micro-optimize QuadEdge and Smooth.
- Make original face edges bold in wireframe
- Textures
- Texture coordinate smoothing
- Animate between levels
- Sprites at fixed screen size

# Edge based wireframe code

```javascript
this.makeWireGeometry = function() {
    var qe = this.qe;
    var edges = qe.edges;
    var edgeCount = edges.length;
    var verts = this.verts;

    var geometry = new THREE.Geometry2(edgeCount * 4);
    var vertices = geometry.vertices;
    var offset = 0;
    var vert;

    for (var i = 0; i < edgeCount; i++) {
        var edge = edges[i];
        vert = verts[edge.vert0];
        vertices[offset++] = vert.x;
        vertices[offset++] = vert.y;
        vertices[offset++] = vert.z;

        vert = verts[edge.vert1];
        vertices[offset++] = vert.x;
        vertices[offset++] = vert.y;
        vertices[offset++] = vert.z;
    }

    return geometry;
}
```

